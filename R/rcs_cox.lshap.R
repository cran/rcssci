#'@title  rcs_cox.lshap
#'
#'@description  restricted cubic splines (RCS) published in SCI.
#'@details Cox models with RCS splines were performed to explore the shape linear or nonlinear(U, inverted U,J,S,L,log,-log,temporary plateau shape)
#'
#'@param data data.frame.Rdata
#'@param knot knot=3-7 or automatic calculate by AIC min
#'@param y outcome=0,1
#'@param time censor time
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
#' rcs_cox.lshap(data=sbpdata, y = "status",x = "sbp",time = "time",
#' prob=0.1,filepath=tempdir())
#'# library(rcssci)
#'# rcs_cox.lshap(knot=4,data=sbpdata, y = "status",x = "sbp",covs=c("age"),
#'# time = "time", prob=0.1,filepath="D:/temp")
#'
#' @export
#' @name rcs_cox.lshap
#'
globalVariables(c('..density..', 'Cairo' ,'aes', 'dplyr' ,'element_blank', 'element_line', 'geom_bar',
                  'geom_density', 'geom_hline' ,'geom_line' ,'geom_point', 'geom_ribbon', 'geom_segment',
                  'geom_text', 'geom_vline', 'ggplot2' ,'ggsave', 'lower', 'patchwork', 'pct', 'plot_layout',
                  'rms', 'scale_x_continuous', 'scale_y_continuous' ,'sec_axis', 'segmented', 'survival',
                  'survminer', 'theme', 'theme_bw', 'upper' ,'yhat','datadist','dd','Surv'))

rcs_cox.lshap<-function(data,knot,y,x,time,covs,prob,filepath,...)
{
  pacman::p_load(rms,ggplot2,survminer,survival,dplyr,segmented,patchwork,Cairo)
  if (!missing(knot)) {warning("please be sure of knot by AIC min(default) or preliminary investigation suggested")}
  if (missing(data)) {stop("data required.")}
  if (missing(x)) {stop("x required.")}
  if (missing(time)) {stop("time required.")}
  if (missing(filepath)) {stop("file path required.")}
  if (missing(prob)) {prob <- 0.5} else {assign("prob",prob)}

  call <- match.call()
  data <- as.data.frame(data)
  y <- y
  x <- x
  time <- time
  if (missing(covs)) {
    covs=NULL
    indf <- dplyr::select(data,y,x,time)
  } else {assign("covs",covs)
    indf <- dplyr::select(data,y,x,time,covs)
  }

  indf[, "y"] <- indf[, y]
  indf[, "x"] <- indf[, x]
  indf[, "time"] <- indf[, time]
  sum(!complete.cases(indf[, c(y, x)]))
  indf <- indf[complete.cases(indf[, c(y, x)]),]
  dd<-NULL
  dd <<- rms::datadist(indf)
  old <- options()
  on.exit(options(old))
  options(datadist = "dd")

  aics <- NULL
  S <- Surv(indf$time,indf$y==1)
  for (i in 3:7) {if (is.null(covs)) {
    formula <- paste0("S~ rcs(x, ", i, ")") }
    else {formula <- paste0("S~ rcs(x, ", i, ")", " + ", paste0(covs, collapse=" + "))
    }
    fit <- rms::cph(as.formula(formula), data=indf, x= TRUE, y= TRUE, tol=1e-25, surv = TRUE)
    summary(fit)
    aics <- c(aics, AIC(fit))
    kn <- seq(3, 7)[which.min(aics)]
  }
  if (missing(knot)){
    knot <- kn
  }
  if (is.null(covs)) {
    formula <- paste0("S~ rcs(x, ", knot, ")")
  }
  else {formula <- paste0("S~ rcs(x, ", knot, ")", " + ", paste0(covs, collapse=" + "))
  }

  model <- rms::cph(as.formula(formula), data=indf, x= TRUE, y= TRUE, tol=1e-25, surv = TRUE)
  model.cox <- model
  phassump <- survival::cox.zph(model, transform="km")
  phresidual <- survminer::ggcoxzph(survival::cox.zph(model,transform="km"))
  anova(model)
  pvalue_all <- anova(model)[1, 3]
  pvalue_nonlin <- round(anova(model)[2, 3],3)
  pre.model <-rms::Predict(model.cox,x,fun=exp,type="predictions",ref.zero=T,conf.int = 0.95,digits=2)

  Q20 <- quantile(indf$x,probs = seq(0,1,0.05))
  segdata <- data.frame(pre.model)
  fit_lm <- lm(yhat~x, data = segdata)
  lshap <- segmented::segmented(fit_lm,seg.Z = ~x, npsi = 1)
  lshap.cutoff <- lshap[["psi"]][1,2]
  dd <<- rms::datadist(indf)
  dd[["limits"]]["Adjust to", "x"] <<- lshap.cutoff
  old <- options()
  on.exit(options(old))
  options(datadist = "dd")

  model <- update(model)
  model.cox <- model
  pre.model <-rms::Predict(model.cox,x,fun=exp,type="predictions",ref.zero=T,conf.int = 0.95,digits=2)


  lshapci <- data.frame(pre.model)
  lshapcicross <- lshapci[which.min(abs(lshapci$x-lshap.cutoff)),]

  ## prob ggplot
  newdf1 <- as.data.frame(dplyr::select(pre.model,x,yhat,lower,upper))
  colnames(newdf1) <- c("x", "y", "lower", "upper")
  min(newdf1[, "x"])
  max(newdf1[, "x"])
  xmin <- min(newdf1[, "x"])
  xmax <- max(newdf1[, "x"])
  min(newdf1[, "lower"])
  max(newdf1[, "upper"])
  ymax1 <- ceiling(max(newdf1[, "upper"]))
  newdf2 <- indf[indf[, "x"] >= xmin & indf[, "x"] <= xmax,]
  breaks <- seq(xmin, xmax, length=20)
  h <- hist(newdf2[, "x"], breaks=breaks, right=TRUE)
  max(h[["counts"]] / sum(h[["counts"]]))
  ymax2 <- 20
  newdf3 <- data.frame(x=h[["mids"]], freq=h[["counts"]], pct=h[["counts"]]/sum(h[["counts"]]))
  freq <- cut(newdf2[, "x"], breaks=breaks, dig.lab=6, right=TRUE)
  as.data.frame(table(freq))
  scale_factor <- ymax2 / ymax1
  xtitle <- x
  ytitle1 <- paste0("HR where the lshap.cutoff for ", x, " is ", sprintf("%.3f", lshap.cutoff))
  ytitle2 <- "Percentage of Population (%)"
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

  ## histogram plot by parameter "lshap.cutoff"
  plot.lshap.type1 <- ggplot2::ggplot() +
    geom_bar(
      data=newdf3, aes(x=x, y=pct*100/scale_factor),
      stat="identity", width=(xmax-xmin)/(length(breaks)-1),
      fill="#f9f7f7", color="grey"
    ) +
    geom_hline(yintercept=1, linetype=2, color="grey") +
    geom_line(data=newdf1, aes(x=x, y=lower), linetype=2, color="#ff9999",size=0.8) +
    geom_line(data=newdf1, aes(x=x, y=upper), linetype=2, color="#ff9999",size=0.8) +
    geom_line(data=newdf1, aes(x=x, y=y), color="#e23e57",size=1) +
    geom_point(aes(x=lshap.cutoff, y=1), color="#e23e57", size=2) +
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
      ytitle1, expand=c(0, 0), limit=c(0, ymax1),
      sec.axis=sec_axis(
        ytitle2, trans=~. * scale_factor,
      )
    ) +
    theme_bw() +
    theme(
      axis.line=element_line(),
      panel.grid=element_blank(),
      panel.border=element_blank()
    )

  ## histogram plot by parameter "lshap.cutoff" with CI shadow
  plot.lshap.type2 <- ggplot2::ggplot() +
    geom_bar(
      data=newdf3, aes(x=x, y=pct*100/scale_factor),
      stat="identity", width=(xmax-xmin)/(length(breaks)-1),
      fill="#f9f7f7", color="grey"
    ) +
    geom_hline(yintercept=1, size=1,linetype=2, color="grey") +
    geom_ribbon(
      data=newdf1, aes(x=x, ymin=lower, ymax=upper),
      fill="#e23e57", alpha=0.1) +
    geom_line(data=newdf1, aes(x=x, y=lower), linetype=0, color="#ff9999",size=0.8) +
    geom_line(data=newdf1, aes(x=x, y=upper), linetype=0, color="#ff9999",size=0.8) +
    geom_line(data=newdf1, aes(x=x, y=y), color="#e23e57",size=1) +
    geom_point(aes(x=lshap.cutoff, y=1), color="#e23e57", size=2) +
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
      ytitle1, expand=c(0, 0), limit=c(0, ymax1),
      sec.axis=sec_axis(
        ytitle2, trans=~. * scale_factor,
      )
    ) +
    theme_bw() +
    theme(
      axis.line=element_line(),
      panel.grid=element_blank(),
      panel.border=element_blank()
    )

  ## simple rcs plot by parameter "lshap.cutoff"
  plot.lshap.type3 <- ggplot2::ggplot() +
    geom_hline(yintercept=1, size=1,linetype=2, color="grey") +
    geom_vline(xintercept=lshap.cutoff,size=1,linetype=1,color = '#d40e8c',alpha=0.3)+
    geom_ribbon(
      data=newdf1, aes(x=x, ymin=lower, ymax=upper),
      fill="#e23e57", alpha=0.1) +
    geom_line(data=newdf1, aes(x=x, y=lower), linetype=0, color="#ff9999",size=0.8) +
    geom_line(data=newdf1, aes(x=x, y=upper), linetype=0, color="#ff9999",size=0.8) +
    geom_line(data=newdf1, aes(x=x, y=y), color="#e23e57",size=1) +
    geom_point(aes(x=lshap.cutoff, y=1), color="#e23e57", size=2) +
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
      ytitle1, expand=c(0, 0), limit=c(0, ymax1)
    ) +
    theme_bw() +
    theme(
      axis.line=element_line(),
      panel.grid=element_blank(),
      panel.border=element_blank()
    )


  ## Density plot by parameter "lshap.cutoff"
  ytitle2 <- "Density"
  ymax2 <- 0.1
  scale_factor2 <- 0.1 / ymax1
  plot.lshap.type4 <- ggplot2::ggplot()+
    geom_density(data = newdf2,
                 mapping =aes(x = x, y = ..density.. /scale_factor2),
                 fill="#f9f7f7", color="grey",
                 alpha =0.8)+
    geom_hline(yintercept=1, size=1,linetype=2, color="grey") +
    geom_ribbon(
      data=newdf1, aes(x=x, ymin=lower, ymax=upper),
      fill="#e23e57", alpha=0.1) +
    geom_line(data=newdf1, aes(x=x, y=lower), linetype=0, color="#ff9999",size=0.8) +
    geom_line(data=newdf1, aes(x=x, y=upper), linetype=0, color="#ff9999",size=0.8) +
    geom_line(data=newdf1, aes(x=x, y=y), color="#e23e57",size=1) +
    geom_point(aes(x=lshap.cutoff, y=1), color="#e23e57", size=2) +
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
      ytitle1, expand=c(0, 0), limit=c(0, ymax1),
      sec.axis=sec_axis(
        ytitle2, trans=~. * scale_factor2,
      )
    ) +
    theme_bw() +
    theme(
      axis.line=element_line(),
      panel.grid=element_blank(),
      panel.border=element_blank()
    )


  fig.lshapall <- plot.lshap.type1+plot.lshap.type2+plot.lshap.type3+plot.lshap.type4+plot_layout(nrow = 2, byrow = TRUE)

  dev.new()
  ggsave("fig.cox_lshapall.PDF", fig.lshapall, width = 14, height =10  , device = cairo_pdf, family = "Times",path=filepath)
  dev.off()

  message.print4 <- list(aics=aics,kn=kn,
                        phassump=phassump,phresidual=phresidual,
                        Q20=Q20,lshapcicross=lshapcicross)
  return(message.print4)
}
















