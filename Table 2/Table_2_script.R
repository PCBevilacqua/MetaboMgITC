devtools::document()
devtools::load_all()

?Kd.app.calc

####Load in Metabolite concentration data####

list.files("Binding_constant_concentration_data")

df = read.csv("Binding_constant_concentration_data/210526_Metabolites.csv")

df$Metabolites

####Calculate apparant Kds for IS = 0.15 and pH = 7.5####

IS = 0.15
pH = 7.5

#Test

Kd.app = Kd.app.calc(df$Metabolites[1], pH, IS)
Kd.app

#Run in a for loop

Kd.app <- c()

for (i in 1:length(df$Metabolites)){
  print(i)
  print(df$Metabolites[i])
  Kd.app[i] = Kd.app.calc(df$Metabolites[i], pH, IS)
}

df$Kd.app = Kd.app

####Print the table and edit by hand####

