read.itc =function(path.to.itc.file,
                   print.peak.integration.graph = FALSE,
                   integration.folder = NA){
  ####Open a connection to a ITC file####
  con = file(path.to.itc.file, "r")

  Temp.Found = FALSE
  reagent.lines = c()
  ####Parses ITC data file into convenient R data frames####
  while ( TRUE ) {
    line = readLines(con, n = 1)
    #print(line)
    if ( length(line) == 0 ) {
      break
    }
    #Finds the temperature on the second line of the ITC file
    if (line != "$ITC"){
      if (Temp.Found == FALSE){
        Temp = as.numeric(as.character(strsplit(line, split = " ")[[1]][2]))
        Temp.Found = TRUE
      }
    }
    #Pulls out injection information
    if (strsplit(line, split = ":")[[1]][1] == "$ADCGainCode"){
      inj.V = c()
      inj.dur = c()
      inj.delay =c()
      while(strsplit(line, split = " ")[[1]][1] != "#"){
        line.vector = strsplit(line, split = " ")[[1]]
        if (length(line.vector) > 1){
          if (line.vector[1] != "#"){
            if (line.vector[1] != "$ADCGainCode:"){
              line.vector = line.vector[which(line.vector != "$")]
              line.vector = line.vector[which(line.vector != ",")]
              inj.V = c(inj.V, as.numeric(as.character(line.vector[1])))
              inj.dur = c(inj.dur, as.numeric(as.character(line.vector[2])))
              inj.delay =c(inj.delay, as.numeric(as.character(line.vector[3])))
              #print(line.vector)
            }
          }
        }
        line = readLines(con, n = 1)
      }
    }
    #Pulls out volumes and reagent concentrations

    if (strsplit(line, split = " ")[[1]][1] == "#"){
      #print(line)
      reagent.lines = c(reagent.lines, strsplit(line, split = " ")[[1]][2])
      #print(reagent.lines)
    }

    #Pulls out delta power data
    if (strsplit(line, split = "")[[1]][1] == "@"){
      inj.n = c()
      Time.sec = c()
      dP = c()
      Temperature = c()
      while(TRUE){
        if (strsplit(line, split = "")[[1]][1] == "@"){
          #print(line)
          inj.N = as.numeric(as.character(gsub("@", "", strsplit(line, split = ",")[[1]][1])))
          #print(inj.N)
          line = readLines(con, n = 1)
        }else{
          #print(line)
          line.vector = strsplit(line, split = ",")[[1]]
          inj.n = c(inj.n, inj.N)
          Time.sec = c(Time.sec, as.numeric(as.character(line.vector[1])))
          dP = c(dP, as.numeric(as.character(line.vector[2])))
          Temperature = c(Temperature, as.numeric(as.character(line.vector[3])))
          line = readLines(con, n = 1)
        }
        if ( length(line) == 0 ) {
          break
        }
      }
    }
  }
  ####Close the connection to the ITC data file####

  close(con)

  ####Parse reagent concentration data####

  V0 = as.numeric(as.character(reagent.lines[4]))
  M0 = as.numeric(as.character(reagent.lines[3]))/1000
  X0 = as.numeric(as.character(reagent.lines[2]))/1000

  ####Make a data frame for the raw differential power####

  df.dP = data.frame(inj.n, Time.sec, dP, Temperature)

  ####make a data frame that summarizes each injection####

  N = 1:length(inj.V)
  df.inj = data.frame(N, inj.V, inj.dur, inj.delay)

  ####Determine the concentration of reagents after each injection####

  V = c()
  X =c()
  M = c()
  dX = c()

  for (i in df.inj$N){
    V[i] = V0 + sum(df.inj$inj.V[1:i])/1000
    X[i] = (sum(df.inj$inj.V[1:i])/1000)*X0/V[i]
    M[i] = V0*M0/V[i]
  }

  df.inj$V = V
  df.inj$X = X
  df.inj$M = M

  ####Determine the heat produced by each injection####

  dQ = c()
  for (i in 1:length(df.inj$N)){
    df <- subset(df.dP, df.dP$inj.n == df.inj$N[i])
    #pracma::polyarea(c(0,1,1,0), c(0, 0, 1, 1))
    #geometry::polyarea(c(0,1,1,0), c(0, 0, 1, 1))
    #geometry::polyarea(df$Time.sec, df$dP)
    dQ[i] = -pracma::polyarea(df$Time.sec, df$dP)
  }

  df.inj$dQ = dQ
  df.inj$dX = (df.inj$inj.V/10^6)*(X0)
  df.inj$dQ.dX = df.inj$dQ/df.inj$dX/(10^9)

  ####Print a graph in order to check peak integration####

  if (print.peak.integration.graph){
    if (length(is.na(integration.folder)) != 1){
      ggplot2::ggplot(df.dP, ggplot2::aes(x = Time.sec, y = dP, color = factor(inj.n), fill = factor(inj.n)))+
        ggplot2::geom_line() +
        ggplot2::geom_polygon() +
        ggplot2::scale_fill_manual(name = "Injection", values = viridis::viridis(length(unique(df.dP$inj.n)))) +
        ggplot2::scale_color_manual(name = "Injection", values = viridis::viridis(length(unique(df.dP$inj.n)))) +
        ggplot2::theme_classic() +
        ggplot2::xlab("Time (sec)") +
        ggplot2::ylab("")
      ggplot2::ggsave(paste(integration.folder, "/", path_to_itc_file,"_integration.png", sep = ""))
    }
  }

  ####Make a list containing the data frames####
  output = list(df.dP, df.inj)
  names(output) = c("dP", "inj")
  print(path.to.itc.file)
  print(names(output)[1])
  print(head(output[[1]]))
  print(names(output)[2])
  print(head(output[[2]]))
  output = output
}


