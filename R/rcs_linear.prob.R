#'@title  rcs_linear.prob
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
#'@examples library(rcssci)
#' rcs_linear.prob(data=sbpdata, y = "sbp",x = "age",
#' prob=0.1,filepath=tempdir())
#'# library(rcssci)
#'# rcs_linear.prob(knot=4,data=sbpdata, y = "sbp",x = "age",
#'# covs=c("gender"),prob=0.1,filepath="D:/temp")
#'
#' @export
#' @name rcs_linear.prob
#'
globalVariables(c('..density..', 'Cairo' ,'aes', 'dplyr' ,'element_blank', 'element_line', 'geom_hline' ,'geom_line' ,'geom_point', 'geom_ribbon', 'geom_segment',
                  'geom_text', 'geom_vline', 'ggplot2' ,'ggsave', 'lower', 'patchwork', 'pct', 'plot_layout',
                  'rms', 'scale_x_continuous', 'scale_y_continuous' ,'sec_axis',  'survival',
                  'survminer', 'theme', 'theme_bw', 'upper' ,'yhat','datadist','dd'))

rcs_linear.prob<-function(data,knot,y,x,covs,prob,filepath,...)
{
  pacman::p_load(rms,ggplot2,survminer,survival,dplyr,patchwork,Cairo)
  if (!missing(knot)) {warning("please be sure of knot by AIC min(default) or preliminary investigation suggested")}
  if (missing(data)) {stop("data required.")}
  if (missing(x)) {stop("x required.")}
  if (missing(prob)) {prob <- 0.5} else {assign("prob",prob)}
  if (missing(filepath)) {stop("file path required.")}

  call <- match.call()
  data <- as.data.frame(data)
  y <- y
  x <- x
  if (missing(covs)) {
    covs=NULL
    indf <- dplyr::select(data,y,x)
  } else {assign("covs",covs)
    indf <- dplyr::select(data,y,x,covs)
  }

  indf[, "y"] <- indf[, y]
  indf[, "x"] <- indf[, x]
  sum(!complete.cases(indf[, c(y, x)]))
  indf <- indf[complete.cases(indf[, c(y, x)]),]
  dd<-NULL
  dd <<- rms::datadist(indf)
  old <- options()
  on.exit(options(old))
  options(datadist = "dd")

  aics <- NULL
  for (i in 3:7) {
    if (is.null(covs)) {formula <- paste0("y~ rcs(x, ", i, ")",paste0(covs, collapse=" + "))
    } else {formula <- paste0("y~ rcs(x, ", i, ")", " + ", paste0(covs, collapse=" + "))}
    fit <- rms::ols(as.formula(formula), data=indf, x=TRUE,tol=1e-25)
    summary(fit)
    aics <- c(aics, AIC(fit))
    kn <- seq(3, 7)[which.min(aics)]
  }
  if (missing(knot)){
    knot <- kn
  }
  if (is.null(covs)) {
    formula <- paste0("y~ rcs(x, ", knot, ")",paste0(covs, collapse=" + "))
  }
  else {
    formula <- paste0("y~ rcs(x, ", knot, ")", " + ", paste0(covs, collapse=" + "))
  }
  model <- rms::ols(as.formula(formula), data=indf, x=TRUE,tol=1e-25)
  model.linear <- model
  anova(model)
  pvalue_all <- anova(model)[1, 5]
  pvalue_nonlin <- round(anova(model)[2, 5],3)
  pre.model <-rms::Predict(model.linear,x,type="predictions",ref.zero=T,conf.int = 0.95,digits=2)

  Q20 <- quantile(indf$x,probs = seq(0,1,0.05))
  probtemp <- prob
  #relatively refvalue = prob
  refvalue <- quantile(indf[, "x"], prob=probtemp)
  dd <<- rms::datadist(indf)
  dd[["limits"]]["Adjust to", "x"] <<- refvalue
  old <- options()
  on.exit(options(old))
  options(datadist = "dd")

  model <- update(model)
  model.linear <- model
  pre.model <-rms::Predict(model.linear,x,type="predictions",ref.zero=T,conf.int = 0.95,digits=2)

  newdf1 <- as.data.frame(dplyr::select(pre.model,x,yhat,lower,upper))
  colnames(newdf1) <- c("x", "y", "lower", "upper")
  min(newdf1[, "x"])
  max(newdf1[, "x"])
  xmin <- min(newdf1[, "x"])
  xmax <- max(newdf1[, "x"])
  min(newdf1[, "lower"])
  max(newdf1[, "upper"])
  ymax1 <- round(max(newdf1[, "upper"]),digits=1)
  ymin1 <- round(min(newdf1[, "lower"]),digits=1)

  xtitle <- x
  ytitle1 <- paste0("continue y where the refvalue for ", x, " is ", sprintf("%.3f", refvalue))

  offsetx1 <- (xmax - xmin) * 0.02
  offsety1 <- ymax1 * 0.02
  labelx1 <- xmin + (xmax - xmin) * 0.2
  labely1 <- ymax1 * 0.9
  label1 <- paste0("Estimation", "\n", "95% CI")
  labelx2 <- xmin + (xmax - xmin) * 0.7
  labely2 <- ymax1 * 0.9
  label2 <- paste0(
    "P-overall = ",
    ifelse(pvalue_all < 0.001, "< 0.001", sprintf("%.3f", pvalue_all)),
    "\nP-non-linear = ",
    ifelse(pvalue_nonlin < 0.001, "< 0.001", sprintf("%.3f", pvalue_nonlin))
  )

  plot.prob.type1 <- ggplot2::ggplot() +
    geom_hline(yintercept=0, linetype=2, color="grey") +
    geom_line(data=newdf1, aes(x=x, y=lower), linetype=2, color="#ff9999",size=0.8) +
    geom_line(data=newdf1, aes(x=x, y=upper), linetype=2, color="#ff9999",size=0.8) +
    geom_line(data=newdf1, aes(x=x, y=y), color="#e23e57",size=1) +
    geom_point(aes(x = refvalue, y=0), color="#e23e57", size=2) +
    geom_segment(
      aes(
        x=c(labelx1-offsetx1*5, labelx1-offsetx1*5),
        xend=c(labelx1-offsetx1, labelx1-offsetx1),
        y=c(labely1+offsety1, labely1-offsety1),
        yend=c(labely1+offsety1, labely1-offsety1)
      ),
      linetype=c(1, 2), color=c("#e23e57", "black")
    ) +
    geom_text(aes(x=labelx1, y=labely1, label=label1), hjust=0) +
    geom_text(aes(x=labelx2, y=labely2, label=label2), hjust=0) +
    scale_x_continuous(xtitle, expand=c(0, 0.01), limit=c(xmin, xmax)) +
    scale_y_continuous(
      ytitle1, expand=c(0, 0), limit=c(ymin1, ymax1),
    ) +
    theme_bw() +
    theme(
      axis.line=element_line(),
      panel.grid=element_blank(),
      panel.border=element_blank()
    )

  plot.prob.type2 <- ggplot2::ggplot() +
    geom_hline(yintercept=0, size=1,linetype=2, color="grey") +
    geom_ribbon(
      data=newdf1, aes(x=x, ymin=lower, ymax=upper),
      fill="#e23e57", alpha=0.1) +
    geom_line(data=newdf1, aes(x=x, y=lower), linetype=0, color="#ff9999",size=0.8) +
    geom_line(data=newdf1, aes(x=x, y=upper), linetype=0, color="#ff9999",size=0.8) +
    geom_line(data=newdf1, aes(x=x, y=y), color="#e23e57",size=1) +
    geom_point(aes(x = refvalue, y=0), color="#e23e57", size=2) +
    geom_segment(
      aes(
        x=c(labelx1-offsetx1*5, labelx1-offsetx1*5),
        xend=c(labelx1-offsetx1, labelx1-offsetx1),
        y=c(labely1+offsety1, labely1-offsety1),
        yend=c(labely1+offsety1, labely1-offsety1)
      ),
      linetype=c(1, 2), color=c("#e23e57", "black")
    ) +
    geom_text(aes(x=labelx1, y=labely1, label=label1), hjust=0) +
    geom_text(aes(x=labelx2, y=labely2, label=label2), hjust=0) +
    scale_x_continuous(xtitle, expand=c(0, 0.01), limit=c(xmin, xmax)) +
    scale_y_continuous(
      ytitle1, expand=c(0, 0), limit=c(ymin1, ymax1),
    ) +
    theme_bw() +
    theme(
      axis.line=element_line(),
      panel.grid=element_blank(),
      panel.border=element_blank()
    )

  plot.prob.type3 <- ggplot2::ggplot() +
    geom_hline(yintercept=0, size=1,linetype=2, color="grey") +
    geom_vline(xintercept=refvalue,size=1,linetype=1,color = '#d40e8c',alpha=0.3)+
    geom_ribbon(
      data=newdf1, aes(x=x, ymin=lower, ymax=upper),
      fill="#e23e57", alpha=0.1) +
    geom_line(data=newdf1, aes(x=x, y=lower), linetype=0, color="#ff9999",size=0.8) +
    geom_line(data=newdf1, aes(x=x, y=upper), linetype=0, color="#ff9999",size=0.8) +
    geom_line(data=newdf1, aes(x=x, y=y), color="#e23e57",size=1) +
    geom_point(aes(x = refvalue, y=0), color="#e23e57", size=2) +
    geom_segment(
      aes(
        x=c(labelx1-offsetx1*5, labelx1-offsetx1*5),
        xend=c(labelx1-offsetx1, labelx1-offsetx1),
        y=c(labely1+offsety1, labely1-offsety1),
        yend=c(labely1+offsety1, labely1-offsety1)
      ),
      linetype=c(1, 2), color=c("#e23e57", "black")
    ) +
    geom_text(aes(x=labelx1, y=labely1, label=label1), hjust=0) +
    geom_text(aes(x=labelx2, y=labely2, label=label2), hjust=0) +
    scale_x_continuous(xtitle, expand=c(0, 0.01), limit=c(xmin, xmax)) +
    scale_y_continuous(
      ytitle1, expand=c(0, 0), limit=c(ymin1, ymax1)
    ) +
    theme_bw() +
    theme(
      axis.line=element_line(),
      panel.grid=element_blank(),
      panel.border=element_blank()
    )


  plot.prob.type4 <- ggplot2::ggplot()+
    geom_hline(yintercept=0, size=1,linetype=2, color="grey") +
    geom_ribbon(
      data=newdf1, aes(x=x, ymin=lower, ymax=upper),
      fill="#e23e57", alpha=0.1) +
    geom_line(data=newdf1, aes(x=x, y=lower), linetype=0, color="#ff9999",size=0.8) +
    geom_line(data=newdf1, aes(x=x, y=upper), linetype=0, color="#ff9999",size=0.8) +
    geom_line(data=newdf1, aes(x=x, y=y), color="#e23e57",size=1) +
    geom_point(aes(x = refvalue, y=0), color="#e23e57", size=2) +
    geom_segment(
      aes(
        x=c(labelx1-offsetx1*5, labelx1-offsetx1*5),
        xend=c(labelx1-offsetx1, labelx1-offsetx1),
        y=c(labely1+offsety1, labely1-offsety1),
        yend=c(labely1+offsety1, labely1-offsety1)
      ),
      linetype=c(1, 2), color=c("#e23e57", "black")
    ) +
    geom_text(aes(x=labelx1, y=labely1, label=label1), hjust=0) +
    geom_text(aes(x=labelx2, y=labely2, label=label2), hjust=0) +
    scale_x_continuous(xtitle, expand=c(0, 0.01), limit=c(xmin, xmax)) +
    scale_y_continuous(
      ytitle1, expand=c(0, 0), limit=c(ymin1, ymax1),
    ) +
    theme_bw() +
    theme(
      axis.line=element_line(),
      panel.grid=element_blank(),
      panel.border=element_blank()
    )
  fig.proball <- plot.prob.type1+plot.prob.type2+plot.prob.type3+plot.prob.type4+plot_layout(nrow = 2, byrow = T)

  dev.new()
  ggsave("fig.linear_proball.pdf", fig.proball, width = 14, height =10 , device = cairo_pdf, family = "Times",path=filepath)
  dev.off()

  message.print <- list(aics=aics,kn=kn,Q20=Q20)
  return(message.print)
}




