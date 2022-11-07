#'@title  rcssci_linear
#'
#'@description  restricted cubic splines (RCS) published in SCI.
#'@details linear models with RCS splines were performed to explore the shape linear or nonlinear(U, inverted U,J,S,L,log,-log,temporary plateau shape)
#'
#'@param data data.frame.Rdata
#'@param knot knot=3-7 or automatic calculate by AIC min
#'@param y outcome=0,1
#'@param covs covariables, univariate analysis  without "covs" command, multivariable analysis  with "covs" command
#'@param prob position parameter,range from 0-1
#'@param x  main exposure and X-axis when plotting
#'@param filepath path of plots output.
#'
#' @importFrom grDevices cairo_pdf dev.new dev.off
#' @importFrom graphics hist
#' @importFrom stats AIC anova as.formula complete.cases lm quantile update
#' @import Cairo
#' @import patchwork
#'
#'@return message.print PH assumption and other message
#'@author Zhiqiang Nie, \email{niezhiqiang@@gdph.org.cn}
#'@examples
#'\donttest{library(rcssci)
#' rcssci_linear(data=sbpdata, y = "sbp",x = "age",
#' prob=0.1,filepath=tempdir())}
#'# library(rcssci)
#'# rcssci_linear(knot=4,data=sbpdata, y = "sbp",x = "age",
#'# covs=c("gender"),prob=0.1,filepath="D:/temp")
#'
#' @export
#' @name rcssci_linear
#'
globalVariables(c('..density..', 'Cairo' ,'aes', 'dplyr' ,'element_blank', 'element_line', 'geom_hline' ,'geom_line' ,'geom_point', 'geom_ribbon', 'geom_segment',
                  'geom_text', 'geom_vline', 'ggplot2' ,'ggsave', 'lower', 'patchwork', 'pct', 'plot_layout',
                  'rms', 'scale_x_continuous', 'scale_y_continuous' ,'sec_axis', 'segmented', 'survival',
                  'survminer', 'theme', 'theme_bw', 'upper' ,'yhat','datadist','dd'))

rcssci_linear<-function(data,knot,y,x,covs,prob,filepath,...)
{
  rcs_linear.prob(data= data, knot=knot,y = y,x = x,covs=covs,prob=, filepath=filepath)
  rcs_linear.ushap(data= data,knot=knot,y = y,x = x,covs=covs,prob=, filepath=filepath)
  rcs_linear.nshap(data= data,knot=knot,y = y,x = x,covs=covs,prob=, filepath=filepath)
  rcs_linear.lshap(data= data,knot=knot,y = y,x = x,covs=covs,prob=, filepath=filepath)
}




