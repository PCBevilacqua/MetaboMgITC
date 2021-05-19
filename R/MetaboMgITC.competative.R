#'Generalized ITC data fitting function for competative ITC
#'
#'Fits an ITC experiment to a user specified thermodynamic model to determine thermodynamic statistics. Raw ".itc" formatted data
#'should be read in using "read.itc". Fits the data to a thermodynamic model then prints a graphical summary of the fit.
#'
#'@param cell Object created by "read.itc" for the ITC run containing the macromolecule.
#'@param competator Objecrt created by "read.itc" for the ITC run containing the macromolecule and the competative ligand.
#'@param blank Object created by "read.itc" for the ITC run containing the ligand tirated into buffer.
#'@param competator.blank Object created by "read.itc" for the ITC run containing the ligand tirated into buffer with the competator.
#'@param Thermodynamic.equation The thermodynamic model that you want to fit the data to. Default = "wiseman.isotherm".
#'@param Fit.start Starting parameters for the non-linear regression. Default = = list(H = 2534, K = 2880).
#'@param Remove.injection Injections you want to remove from subsequent analysis. Default = 1 is the standard for ITC experiments. More than one injection can be supplied in a vector.
#'@param Saturation.threshold For low c value ITC experiments fit to Wiseman isotherms, macomolecule saturation ranging from 70% to 90% produce the same answers. Thus, sometimes it is desirable to standardize the saturation threshold between experiments to minimize degredation of a ligand or macromolecule (for example with ATP binding Mg). Default = FALSE. If set to 0.8, this function will only fit data required to reach 80% macromolecule saturation.
#'@return A list containing summary statistics for the fit and a graphical depiction of the fits, the raw ITC curve, and the summary statistics.
#' @export
MetaboMgITC.compete = function(cell,
                               competator,
                               blank,
                               blank.competator,
                               Thermodynamic.equation = "Wiseman.isotherm",
                               Fit.start = list(V0 = 1.4247, H = 2534, K = 2880, Hc = 2534, Kc = 2880, c = 0, h = 0, Mg.contaminate = 1),
                               Remove.injection = 1,
                               Saturation.threshold = FALSE,
                               Save.path.and.prefix = FALSE){
  ####Subtract out the background####

  df1 = cell$inj
  df2 = competator$inj
  df3 = blank.competator$inj
  df.blank = blank$inj
  df1$dQ.dX = df1$dQ.dX - df.blank$dQ.dX
  df2$dQ.dX = df2$dQ.dX - df.blank$dQ.dX
  df3$dQ.dX = df3$dQ.dX - df.blank$dQ.dX
  df1$Sample = "No competator"
  df2$Sample = "Competator"
  df3$Sample = "Competator no Mg"
  df1$N.Mg.contaminate = 0
  df2$N.Mg.contaminate = 1
  df3$N.Mg.contaminate = 1

  ####Define the function####

  if (Thermodynamic.equation == "Wiseman.isotherm"){
    Thermo.eq =  function(V, X, M, C, N.Mg.contaminate, V0, H = 3.1, K = 23328.7, Hc = 0, Kc = 900, c = 0, h = 0, Mg.contaminate = 0.0005){
      Mcontaminate = (Mg.contaminate*N.Mg.contaminate*V0)/(V)
      x = X/(M + Mcontaminate)
      Kapp = K*(1 + c*Kc*C)/(1 + Kc*C)
      Happ = H - ((Hc*Kc*C)/(1+Kc*C)) + ((Hc + h)*(x*Kc*C)/(1 + c*Kc*C))
      r = 1/(Kapp*(M + Mcontaminate))
      y = 0.5 + ((1 - x - r)/(2*sqrt(((1 + x + r)^2)-(4*x))))
      dQ.dX = Happ*V*y
    }
  }

  ####Fit the no competator data####

  df = dplyr::bind_rows(df1)

  df = df[-which(df$N == Remove.injection),]

  df$V0 = 1.4247

  fit1 = nls(dQ.dX ~ Thermo.eq(V, X, M, C, N.Mg.contaminate, V0, H, K, Hc = 0, Kc = 0, c = 0, h = 0, Mg.contaminate = 0), df,
            start = list(H = 3, K = 23000), trace = TRUE, control = nls.control(warnOnly = TRUE),
            algorithm = "port", lower = 0)

  ####Fit the competator blank data####

  df = dplyr::bind_rows(df3)

  df = df[-which(df$N == Remove.injection),]

  df$V0 = 1.4247
  df$H = coef(fit1)[1]
  df$K = coef(fit1)[2]

  fit = nls(dQ.dX ~ Thermo.eq(V, X, M, C, N.Mg.contaminate, V0, H, K, Hc, Kc , c, h, Mg.contaminate), df,
            start = list( Hc = 0.01, Kc = 1000, c = 0, h = 0, Mg.contaminate = 0.0006), trace = TRUE, control = nls.control(warnOnly = TRUE),
            algorithm = "port", lower = 0)

  ####Fit the competator Mg data and the competator blank data####

  ####Fit all the data####

  df = dplyr::bind_rows(df1, df2, df3)

  df = df[-which(df$N == Remove.injection),]

  df$V0 = 1.4247

  fit = nls(dQ.dX ~ Thermo.eq(V, X, M, C, N.Mg.contaminate, V0, H, K, Hc, Kc , c, h, Mg.contaminate),
            df,
            start = list(H = coef(fit1)[1],
                         K = coef(fit1)[2],
                         Hc = coef(fit)[1],
                         Kc = coef(fit)[2],
                         c = coef(fit)[3],
                         h = coef(fit)[4],
                         Mg.contaminate = coef(fit)[5]),
            trace = TRUE,
            control = nls.control(warnOnly = TRUE),
            algorithm = "port", lower = 0)

  coef(fit)

  M = df$M + coef(fit)[7]*df$N.Mg.contaminate*1.4247/df$V

  df$fit = predict(fit)

  df$M = M

  ggplot(df, aes(x = X/M, y = dQ.dX, color = Sample)) + geom_point() + geom_line(mapping = aes(y = fit))

  ####Estimate the Mg contaminate#####

  head(df3)

  df3 = df3[-1,]

  y = c()
  Mcontaminate = c()
  x = c()

  for (i in 1:length(df3$N)){
    Mcontaminate[i] = 0.0005*1.4247/df3$V[i]
    x[i] = df3$X[i]/Mcontaminate[i]
    y[i] = Thermo.eq(df3$V[i], df3$X[i], df3$M[i], C = df3$C[i], N.Mg.contaminate = 1, V0 = 1.4247, H = 3.100, K = 23328.7, Hc = 0.01, Kc = 900, c = 0, h = 0, Mg.contaminate = 0.0005)
  }

  df3$M = Mcontaminate
  df3$y = y

  ggplot(df3, aes(x = X/M)) +
    geom_point(mapping = aes(y = dQ.dX)) +
    geom_line(mapping = aes(y = y))


  df3$M = 0

  df3$V0 = 1.4247
  df3$H = 3.06
  df3$K = 23328.7

  head(df3)

  fit = nls(dQ.dX ~ Thermo.eq(V, X, M, C, N.Mg.contaminate, V0, H, K, Hc, Kc , c, h, Mg.contaminate), df3,
            start = list( Hc = 0.01, Kc = 1000, c = 0, h = 0, Mg.contaminate = 0.0006), trace = TRUE, control = nls.control(warnOnly = TRUE),
            algorithm = "port", lower = 0)

  predict(fit)

  coef(fit)

  df = dplyr::bind_rows(df2, df3)

  df = df[-which(df$N == Remove.injection),]

  df$V0 = 1.4247
  df$H = 3.06
  df$K = 23328.7

  head(df)

  fit = nls(dQ.dX ~ Thermo.eq(V, X, M, C, N.Mg.contaminate, V0, H, K, Hc, Kc , c, h, Mg.contaminate), df3,
            start = list(Hc = coef(fit)[1], Kc = coef(fit)[2], c = coef(fit)[3], h = coef(fit)[4], Mg.contaminate = coef(fit)[5]), trace = TRUE, control = nls.control(warnOnly = TRUE),
            algorithm = "port", lower = 0)



  coef(fit)

  df = dplyr::bind_rows(df1, df2, df3)

  head(df)

  df = df[-which(df$N == Remove.injection),]

  df$V0 = 1.4247

  fit = nls(dQ.dX ~ Thermo.eq(V, X, M, C, N.Mg.contaminate, V0, H, K, Hc, Kc , c, h, Mg.contaminate), df3,
            start = list(H = 3.06, K = 23328.7, Hc = coef(fit)[1], Kc = coef(fit)[2], c = coef(fit)[3], h = coef(fit)[4], Mg.contaminate = coef(fit)[5]), trace = TRUE, control = nls.control(warnOnly = TRUE),
            algorithm = "port", lower = 0)

  coef(fit)






  ####Remove bad injections####



  ####Fit data####

  #?nls

  y = c()

  for (i in 1:length(df$N)){
    y[i] = Thermo.eq(df$V[i], df$X[i], df$M[i], df$C[i], df$N.Mg.contaminate[i], V0 = 1.4247, H = 2534, K = 2880, Hc = 2534, Kc = 2880, c = 0, h = 0, Mg.contaminate = 0.0001)
  }

  plot(df$X/df$M, y)

  fit = nls(dQ.dX ~ Thermo.eq(V, X, M, C, N.Mg.contaminate, V0 = 1.4247, H, K, Hc, Kc, c, h, Mg.contaminate), df, start = Fit.start, trace = TRUE, control = nls.control(warnOnly = TRUE),
            algorithm = "port", lower = 0)

  ####Make a heat plot####

  fit.predict = predict(fit)
  fit.predict = c(fit.predict, rep(NA, length(df$N) - length(fit.predict)))

  df$fit = fit.predict

  Heat.plot <- ggplot2::ggplot(df, ggplot2::aes(x = X/M, y = dQ.dX, color = Sample)) +
    ggplot2::geom_point() +
    ggplot2::geom_line(mapping = ggplot2::aes(x = X/M, y = fit)) +
    ggplot2::theme_classic() +
    ggplot2::scale_color_manual(values = c("red", "black")) +
    ggplot2::geom_vline(xintercept = 1) +
    ggplot2::xlab("X/M") +
    ggplot2::ylab("kcal/mol of X") +
    ggplot2::theme(legend.position = "none")

  ####Make a differential power plot####

  cell$dP$Sample = "No competator"
  competator$dP$Sample = "Competator"
  df.dPower = dplyr::bind_rows(cell$dP, competator$dP)

  dP.plot <- ggplot2::ggplot(df.dPower, ggplot2::aes(x = Time.sec/60, y = dP, color = Sample)) +
    ggplot2::geom_line() +
    ggplot2::theme_classic() +
    ggplot2::xlab("Time (min)") +
    ggplot2::facet_wrap(~Sample, ncol = 1) +
    ggplot2::scale_color_manual(values = c("red", "black")) +
    ggplot2::ylab(paste("\U003BC", "cal/sec", sep = "")) +
    ggplot2::theme(legend.position = "none")

  ####Tabulate the fit results####

  model = Thermodynamic.equation
  n = 1
  Temperature = mean(cell$dP$Temperature)
  K = coef(fit)[2]
  dG = Gibbs.free.energy(K, Temperature)
  dH = coef(fit)[1]
  dS = Entropy(dG, dH, Temperature)
  Kc = coef(fit)[4]
  dGc = Gibbs.free.energy(Kc, Temperature)
  dHc = coef(fit)[3]
  dSc = Entropy(dGc, dHc, Temperature)

  Sample = c(rep("EDTA", 6), rep("Competator", 6))
  Parameter = c("n", "Temp.", "K", "dG", "dH", "dS",
                "n", "Temp.", "K", "dG", "dH", "dS")
  Number = c(n, Temperature, K, dG, dH, dS,
             n, Temperature, Kc, dGc, dHc, dSc)

  df.table = data.frame(Sample, Parameter, round(Number, digits = 2))
  df.table$Parameter = factor(df.table$Parameter, levels = c("n", "Temp.", "K", "dG", "dH", "dS"))
  colnames(df.table)[3] = "Number"

  results.table = ggplot2::ggplot(df.table, ggplot2::aes(x = Parameter, y = Sample, label = round(Number, digits = 2))) +
    ggplot2::geom_tile(fill = "white", color = "black") +
    ggplot2::geom_text() +
    ggplot2::theme_classic() +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   axis.line = ggplot2::element_blank())

  ####Make a summary graph####

  final.plot = cowplot::plot_grid(dP.plot, Heat.plot, results.table, ncol = 1, rel_heights = c(2.5,1.5,1))

  ####Save the data####

  if (Save.path.and.prefix != FALSE){
    ggplot2::ggsave(paste(Save.path.and.prefix, "_results_plot.png", sep = ""), width = 6, height = 10)
    write.csv(df.table, paste(Save.path.and.prefix, "_results_table.csv", sep = ""), row.names = FALSE)
  }

  ####Output####
  print(df.table)
  output = list(df.table, final.plot)
  names(output) = c("Table", "Plot")
  output = output
}
