#' @title ggplot theme for minimalist plotting
#'
#' @description Theme used for plotting SQUID component plots. Modification of theme_classic from original ggplot2 package, with consistent font sizing across text elements and with black text, rather than dark grey.
#'
#' @param base_size base font size, given in pts
#'
#' @return An object of class \code{\link[ggplot2]{theme}()}.
#'
#' @examples
#' library("ggplot2")
#'
#' ggplot(mtcars) +
#'   aes(x = wt, y = mpg, colour = factor(gear)) +
#'   geom_point() +
#'   theme_squid()
#' @export

theme_squid <- function(base_size = 12) {
  theme_classic() +
    theme(axis.text = element_text(color = "black", size = base_size),
          axis.ticks = element_line(color = "black"),
          plot.title = element_text(size = base_size),
          legend.text = element_text(size = base_size),
          strip.text = element_text(size = base_size),
          strip.background = element_rect(fill = "grey90", color = "black", linewidth = 0.5))
}
