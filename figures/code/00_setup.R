suppressPackageStartupMessages({
  library(tidyverse)
  library(ggpubr)
  # Point labels that shift aside instead of overplotting each other (Figure 4b).
  library(ggrepel)
})

find_root <- function() {
  root <- tryCatch(
    system("git rev-parse --show-toplevel", intern = TRUE, ignore.stderr = TRUE),
    error = function(e) character(0)
  )
  if (length(root) == 1 && dir.exists(root)) root else getwd()
}

ROOT     <- find_root()
FIG_ROOT <- file.path(ROOT, "figures")
DIR_CODE <- file.path(FIG_ROOT, "code")
DIR_FIGS <- file.path(FIG_ROOT, "figs")
DIR_TABS <- file.path(FIG_ROOT, "tabs")

walk(c(DIR_FIGS, DIR_TABS), dir.create, showWarnings = FALSE, recursive = TRUE)

# ggarrange (via cowplot) measures text on a null PDF device; match it to ragg output fonts
cowplot::set_null_device("agg")

run_task <- function(run) {
  run |>
    str_remove("^\\d{4}-\\d{2}-\\d{2}_") |>
    str_remove("_(xgb|wt_wt).*$") |>
    str_replace("^permeability_", "perm. ") |>
    str_replace("^nf$", "non-fouling") |>
    str_replace_all("_", " ")
}

# Paul Tol "muted" (colour-blind safe); one fixed colour per property across all figures
PAL_TASK <- c(
  "binding affinity" = "#332288",
  "half life"        = "#88CCEE",
  "hemolysis"        = "#CC6677",
  "non-fouling"      = "#117733",
  "perm. caco2"      = "#DDCC77",
  "perm. pampa"      = "#999933",
  "perm. penetrance" = "#44AA99",
  "solubility"       = "#AA4499",
  "toxicity"         = "#882255"
)

# MuriVerse evidence-base colours, drawn from PAL_TASK so the overview, Figure 2 and the
# supplementary figures share a palette. Defined once here because an assay setting
# must mean the same colour in every figure that draws it.
# Mixing toward white lifts lightness and drops saturation while holding hue exactly, so the
# pastel palette below is the same three colours, not three new ones.
# The MuriVerse figures (2, 3 and 4) share one blue/green/red triple, defined once here so the
# three cannot drift apart: Figures 3 and 4 give it to the three model branches, Figure 2 to
# the evidence base those branches are trained on. These are Paul Tol muted colours, and
COL_MV_BLUE  <- "#88CCEE"
COL_MV_GREEN <- "#117733"
COL_MV_RED   <- "#CC6677"

# Figure 3e grades how much of each molecule its SMILES covers, so it wants one ramp rather
# than three unrelated hues: full green for the whole molecule, teal for the peptide alone,
# grey for nothing. Teal rather than the sand #DDCC77, which 3f spends on the validation
# split immediately to its right.
COL_MV_TEAL <- "#44AA99"

# The rodent species palette. Figure 2 describes the cohort and Figure 3 models it, so both
# split by species and the two pastel hues are defined once here rather than in either
# figure: mouse is always this pastel navy and rat always this pastel forest green.
MV_MOUSE <- "#4169E1"   # royal blue (fig2, fig3 species palette)
MV_RAT   <- "#2E8B57"   # medium sea green

# Figure 2 spends the blue and the red on large filled areas rather than on the lines and
# points Figures 3 and 4 use them for, so it takes them a shade darker and then draws them
# partly transparent: the darker base keeps the hue from washing out at scale, the
# transparency keeps the result pastel. The green is used as given.
MV_DARKEN <- 0.20     # blue: pulled down from the Figure 3/4 tone
MV_BLOOD  <- 0.70     # red: pulled this far toward blood red, away from rose
MV_PASTEL <- 0.30     # green: lifted this far toward white
MV_FILL_ALPHA <- 0.75

mix_col <- function(a, b, t) {
  m <- grDevices::col2rgb(a) / 255 * (1 - t) + grDevices::col2rgb(b) / 255 * t
  grDevices::rgb(m["red", ], m["green", ], m["blue", ])
}
darken <- function(colour, amount = MV_DARKEN) mix_col(colour, "#000000", amount)

COL_IN_VIVO  <- darken(COL_MV_BLUE)                            # -> #6DA3BE, drawn #92BACE
COL_IN_VITRO <- mix_col(COL_MV_GREEN, "#FFFFFF", MV_PASTEL)    # -> #58A070, drawn #82B894
COL_ANALYSIS <- mix_col(darken(COL_MV_RED), "#8A0303", MV_BLOOD)  # -> #921B1F, drawn #AD5457
# ex vivo sits outside the triple: every remaining Tol muted colour is already spoken for, and
# this salmon is far enough from COL_MV_RED to stay separable where the two appear together
COL_EX_VIVO  <- "#F39B7F"

# The red fill (COL_ANALYSIS as drawn, at MV_FILL_ALPHA) is defined for the audit panels
# in fig2.R (source concentration, homology clusters) which are kept available for ad-hoc
# use. Figure 2's funnel and filtering panels are coloured by species (royal blue mouse,
# green rat); Supplementary Figure 1 uses it to mark the kept side of each gate.
MV_RED_LIGHT <- "#F6D5D3"
MV_RED_DARK  <- mix_col(COL_ANALYSIS, "#FFFFFF", 1 - MV_FILL_ALPHA)

PAL_SETTING <- c(`in vivo` = COL_IN_VIVO, `ex vivo` = COL_EX_VIVO,
                 `in vitro` = COL_IN_VITRO, other = "grey60", unknown = "grey60")

# Fades used for ordered categories, scaled so the solid end matches the MV_FILL_ALPHA the
# rest of Figure 2 is drawn at. The 0.4 floor is set by COL_ANALYSIS being a mid-tone rose
# rather than a dark red: below it the palest step stops reading on white.
FADE_RANGE <- c(0.4, 1) * MV_FILL_ALPHA
fade_pal <- function(colour, n, range = FADE_RANGE) {
  scales::alpha(colour, seq(range[1], range[2], length.out = n))
}

# Binary-classification helpers for the half-life ranking sensitivity analysis.
roc_curve <- function(y, p) {
  o <- order(p, decreasing = TRUE); y <- y[o]
  tibble(fpr = c(0, cumsum(y == 0) / sum(y == 0)),
         tpr = c(0, cumsum(y == 1) / sum(y == 1)))
}

pr_curve <- function(y, p) {
  o <- order(p, decreasing = TRUE); y <- y[o]
  tp <- cumsum(y == 1)
  tibble(recall = tp / sum(y == 1), precision = tp / seq_along(y))
}

auc_trapz <- function(x, y) sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)

# Half-life axes are drawn in clock time (1 m ... 480 h), never in log10 units. HL_BREAKS_H is the
# menu of ticks; each panel picks the few that fit. hl_lab() writes 0.25 as "15 m" and 6 as "6 h".
HL_BREAKS_H <- c(1 / 60, 5 / 60, 0.25, 0.5, 1, 2, 4, 6, 12, 24, 48, 72, 96, 120, 240, 480)
hl_lab <- function(h) ifelse(h < 1, paste0(round(h * 60), " m"), paste0(h, " h"))
# for an axis that holds log10(hours) as a plain number
scale_x_loghl <- function(breaks = c(0.25, 1, 6, 24), ...)
  scale_x_continuous(breaks = log10(breaks), labels = hl_lab(breaks), ...)
scale_y_loghl <- function(breaks = c(0.25, 1, 6, 24), ...)
  scale_y_continuous(breaks = log10(breaks), labels = hl_lab(breaks), ...)
# for an axis that holds hours on a log10 scale
scale_x_hl <- function(breaks = c(0.25, 1, 6, 24), ...)
  scale_x_log10(breaks = breaks, labels = hl_lab(breaks), ...)
scale_y_hl <- function(breaks = c(0.25, 1, 6, 24), ...)
  scale_y_log10(breaks = breaks, labels = hl_lab(breaks), ...)

# Bottom-axis tick labels always run down and to the right at 45 degrees: a negative angle
# descends left-to-right, and hjust = 0 with vjust = 1 anchors each label's start at its tick.
theme_set(
  theme_pubr(base_size = 10) +
    theme(legend.position = "right", strip.background = element_blank(),
          plot.margin = margin(14, 16, 4, 8),
          axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1))
)

save_fig <- function(plot, name, width = 180, height = 160) {
  ggsave(file.path(DIR_FIGS, paste0(name, ".pdf")), plot,
         width = width, height = height, units = "mm", device = cairo_pdf, bg = "white")
  ggsave(file.path(DIR_FIGS, paste0(name, ".png")), plot,
         width = width, height = height, units = "mm", dpi = 300, device = ragg::agg_png,
         bg = "white")
  invisible(plot)
}
