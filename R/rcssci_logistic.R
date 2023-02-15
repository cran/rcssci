#'@title  rcssci_logistic
#'
#'@description  restricted cubic splines (RCS) published in SCI.
#'@details logistic models with RCS splines were performed to explore the shape linear or nonlinear(U, inverted U,J,S,L,log,-log,temporary plateau shape)
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
#' rcssci_logistic(data=sbpdata, y = "status",x = "sbp",
#' prob=0.1,filepath=tempdir())}
#'# library(rcssci)
#'# rcssci_logistic(knot=4,data=sbpdata, y = "status",x = "sbp",
#'# covs=c("age","gender"),prob=0.1,filepath="D:/temp")
#'
#' @export
#' @name rcssci_logistic
#'
globalVariables(c('..density..', 'Cairo' ,'aes', 'dplyr' ,'element_blank', 'element_line', 'geom_bar',
                  'geom_density', 'geom_hline' ,'geom_line' ,'geom_point', 'geom_ribbon', 'geom_segment',
                  'geom_text', 'geom_vline', 'ggplot2' ,'ggsave', 'lower', 'patchwork', 'pct', 'plot_layout',
                  'rms', 'scale_x_continuous', 'scale_y_continuous' ,'sec_axis', 'segmented', 'survival',
                  'survminer', 'theme', 'theme_bw', 'upper' ,'yhat','datadist','dd'))

rcssci_logistic<-function(data,knot,y,x,covs,prob,filepath,...)
{
  if (!missing(knot)) {warning("please be sure of knot by AIC min(default) or preliminary investigation suggested")}
  rcs_logistic.prob(data= data, knot=knot,y = y,x = x,covs=covs,prob=prob, filepath=filepath)
  rcs_logistic.ushap(data= data,knot=knot,y = y,x = x,covs=covs,prob=prob, filepath=filepath)
  rcs_logistic.nshap(data= data,knot=knot,y = y,x = x,covs=covs,prob=prob, filepath=filepath)
  rcs_logistic.lshap(data= data,knot=knot,y = y,x = x,covs=covs,prob=prob, filepath=filepath)
}




