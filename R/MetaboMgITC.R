MetaboMgITC = function(cell,
                       blank,
                       Thermodynamic.equation = "Wiseman.isotherm",
                       Fit.start = list(H = 2534, K = 2880),
                       Remove.injection = 1,
                       Saturation.threshold = FALSE,
                       Save.path.and.prefix = FALSE){
  ####Subtract out the background####

  df = cell$inj
  df.blank = blank$inj
  df$dQ.dX = df$dQ.dX - df.blank$dQ.dX

  ####Define the function####

  if (Thermodynamic.equation == "Wiseman.isotherm"){
    Thermo.eq = Wiseman.isotherm
  }

  ####Remove bad injections####

  df = df[-Remove.injection,]

  ####Fit data####

  fit = nls(dQ.dX ~ Thermo.eq(V, X, M, H, K), df, start = Fit.start)

  ####Calculate saturation####

  K = coef(fit)[2]
  Xf = df[length(df$N),]$X[1]
  Mf = df[length(df$N),]$M[1]
  a = K
  b = K*Xf - K*Mf + 1
  c = -Mf
  M.free = (-b + sqrt(((b^2)-(4*a*c))))/(2*a)
  Sat = 100 - (100*M.free/Mf)

  ####Recursively fit until you are only fitting data points to get to the saturation threshold####

  if (Saturation.threshold != FALSE){
    Remove.saturated.points = TRUE
    Started.recursion =FALSE
    while (Remove.saturated.points){
      if (Started.recursion){
        Fit.start = list(H = coef(fit)[1], K = coef(fit)[2])
        fit = nls(dQ.dX ~ Thermo.eq(V, X, M, H, K), df.trimed, start = Fit.start)
      }else{
        df.trimed = df
      }
      #Calculate final saturation
      K = coef(fit)[2]
      Xf = df.trimed[length(df.trimed$N),]$X[1]
      Mf = df.trimed[length(df.trimed$N),]$M[1]
      a = K
      b = K*Xf - K*Mf + 1
      c = -Mf
      M.free = (-b + sqrt(((b^2)-(4*a*c))))/(2*a)
      Sat = 100 - (100*M.free/Mf)
      #Calculate saturation at each injection
      Sat.n = c()
      for (i in 1:length(df$N)){
        Xf = df$X[i]
        Mf = df$M[i]
        a = K
        b = K*Xf - K*Mf + 1
        c = -Mf
        M.free = (-b + sqrt(((b^2)-(4*a*c))))/(2*a)
        Sat.n[i] = 100 - (100*M.free/Mf)
      }
      sat.injections = which(Sat.n/100 > Saturation.threshold)
      if (length(df$N) - length(sat.injections) == length(df.trimed$N)){
        Remove.saturated.points =  FALSE
      }else{
        df.trimed = df[-sat.injections,]
        Started.recursion = TRUE
      }
    }
  }

  ####Make a heat plot####

  fit.predict = predict(fit)
  fit.predict = c(fit.predict, rep(NA, length(df$N) - length(fit.predict)))

  df$fit = fit.predict

  Heat.plot <- ggplot2::ggplot(df, ggplot2::aes(x = X/M, y = dQ.dX)) +
    ggplot2::geom_point() +
    ggplot2::geom_line(mapping = ggplot2::aes(x = X/M, y = fit)) +
    ggplot2::theme_classic() +
    ggplot2::xlab("X/M") +
    ggplot2::ylab("kcal/mol of X")

  ####Make a differential power plot####

  dP.plot <- ggplot2::ggplot(cell$dP, ggplot2::aes(x = Time.sec/60, y = dP)) +
    ggplot2::geom_line() +
    ggplot2::theme_classic() +
    ggplot2::xlab("Time (min)") +
    ggplot2::ylab(paste("\U003BC", "cal/sec", sep = ""))

  ####Tabulate the fit results####

  model = Thermodynamic.equation
  n = 1
  Temperature = mean(cell$dP$Temperature)
  K = coef(fit)[2]
  dG = Gibbs.free.energy(K, Temperature)
  dH = coef(fit)[1]
  dS = Entropy(dG, dH, Temperature)


  Parameter = c("n", "Temp.", "K", "Saturation", "dG", "dH", "dS")
  Number = c(n, Temperature, K, Sat, dG, dH, dS)

  df.table = data.frame(Parameter, round(Number, digits = 2))
  df.table$Parameter = factor(df.table$Parameter, levels = Parameter)
  colnames(df.table)[2] = "Number"

  results.table = ggplot2::ggplot(df.table, ggplot2::aes(x = Parameter, y = 1, label = round(Number, digits = 2))) +
    ggplot2::geom_tile(fill = "white", color = "black") +
    ggplot2::geom_text() +
    ggplot2::theme_classic() +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::theme(axis.ticks.x = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_blank(),
          axis.line = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_blank())

  ####Make a summary graph####

  final.plot = cowplot::plot_grid(dP.plot, Heat.plot, results.table, ncol = 1, rel_heights = c(2.5,2.5,1))

  ####Save the data####

  if (Save.path.and.prefix != FALSE){
    ggplot2::ggsave(paste(Save.path.and.prefix, "results_plot.png", sep = ""), width = 6, height = 10)
    write.csv(df.table, paste(Save.path.and.prefix, "results_table.csv", sep = ""), row.names = FALSE)
  }

  ####Output####
  print(df.table)
  output = list(df.table, final.plot)
  names(output) = c("Table", "Plot")
  output = output
}
