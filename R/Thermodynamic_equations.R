#'Function that describes heat produced by each injection in an ITC experiment
#'
#'Function that describes heat produced by each injection in an ITC experiment when there is 1 to 1 binding of ligand to a macromolecule and a competator.
#'
#'@param V Cell volume in Liters
#'@param X Ligand concentration in a cell
#'@param M Macromolecule in a cell
#'@param H Enthalpy of a ligand binding to a macromolecule
#'@param K Molar affinity constant of a ligand binding to a macromolecule
#'@return Heat produced by an injection
#' @export
Wiseman.isotherm.competative = function(V, X, M, C, H = 2534, K = 2880, Hc = 2534, Kc = 2880, c = 0, h = 0){
  x = X/M
  Kapp = K*(1 + c*Kc*C)/(1 + Kc*C)
  Happ = H - ((Hc*Kc*C)/(1+Kc*C)) + ((Hc + h)*(x*Kc*C)/(1 + c*Kc*C))
  r = 1/(Kapp*M)
  y = 0.5 + ((1 - x - r)/(2*sqrt(((1 + x + r)^2)-(4*x))))
  dQ.dX = Happ*V*y
}

#'Function that describes heat produced by each injection in an ITC experiment
#'
#'Function that describes heat produced by each injection in an ITC experiment when there is 1 to 1 binding of ligand to a macromolecule.
#'
#'@param V Cell volume in Liters
#'@param X Ligand concentration in a cell
#'@param M Macromolecule in a cell
#'@param H Enthalpy of a ligand binding to a macromolecule
#'@param K Molar affinity constant of a ligand binding to a macromolecule
#'@return Heat produced by an injection
#' @export
Wiseman.isotherm = function(V, X, M, H = 2534, K = 2880){
  x = X/M
  r = 1/(K*M)
  y = 0.5 + ((1 - x - r)/(2*sqrt(((1 + x + r)^2)-(4*x))))
  dQ.dX = H*V*y
}

#'Function that calculates the Gibbs free energy of a chemical reaction
#'
#'@param K Molar affinity constant of a ligand binding to a macromolecule
#'@param Temperature Temperature in celcius
#'@param R Gas constant. Default = 0.0019872 kcal/mol/K
#'@return Gibbs free energy of a chemical reaction
#' @export
Gibbs.free.energy = function(K, Temperature, R = 0.0019872){
  output = -R*(Temperature + 273.15)*log(K)
}

#'Function that calculates the Entropy of a chemical reaction
#'
#'@param dG Gibbs free energy of a chemical reaction
#'@param dH Enthalpy of a chemical reaction
#'@param Temperature Temperature in celcius
#'@return Entropy of a chemical reaction
#' @export
Entropy = function(dG, dH, Temperature){
  output = (dH - dG)/Temperature
}

#'Function that calculates the saturation of a macromolecule in an ITC experiment
#'
#'Function that describes heat produced by each injection in an ITC experiment when there is 1 to 1 binding of ligand to a macromolecule.
#'
#'@param K Molar affinity constant of a ligand binding to a macromolecule
#'@param M Macromolecule in a cell
#'@param X Ligand concentration in a cell
#'@param V Cell volume in Liters
#'@param V0 Original cell volume
#'@return Percent saturation of the macromolecule by bound ligand
#' @export
saturation = function(K, M, X, V, V0){
  M = (V*M)/(1.4 + 0.282)
  X = (0.282*X)/(1.4 + 0.282)
  a = K
  b = K*X - K*M + 1
  c = -M
  M.free = (-b + sqrt(((b^2)-(4*a*c))))/(2*a)
  output = 100 - (100*M.free/M)
  print(output)
}

#'Function that calculates the log10 of a metal ion affinity constant for Ionic strength = 0 using the Debye-Huckel equation
#'
#'@param log10K Molar affinity constant of a metabolite binding to a metal ion
#'@param I ionic strength
#'@param M.charge charge of the metal ion
#'@param X.charge charge of the metabolite
#'@param MX.charge Charge of the metal ion-metabolite complex
#'@param A Debye-Huckel constant for the given temperature. Default A = 0.524 for 25C.
#'@return the log10K at an ionic strength of zero
#' @export
log10K.0.calculator = function(log10K, I, M.charge, X.charge, XM.charge, A = 0.524){
  y.XM = 10^(0.1*(XM.charge^2)*I - (A*(XM.charge^2)*sqrt(I))/(1+sqrt(I)))
  y.M = 10^(0.1*(M.charge^2)*I - (A*(M.charge^2)*sqrt(I))/(1+sqrt(I)))
  y.X = 10^(0.1*(X.charge^2)*I - (A*(X.charge^2)*sqrt(I))/(1+sqrt(I)))
  log10K0 = log10K - log10(y.XM/(y.M*y.X))
}

#'Function that calculates the log10 of a metal ion affinity constant for any Ionic strength  using the Debye-Huckel equation
#'
#'@param log10K.ref Molar affinity constant of a metabolite binding to a metal ion that you are referencing
#'@param I.ref Ionic strength of the reference constant
#'@param I ionic strength for the condition you want to calculate
#'@param M.charge charge of the metal ion
#'@param X.charge charge of the metabolite
#'@param XM.charge Charge of the metal ion-metabolite complex
#'@param A Debye-Huckel constant for the given temperature. Default A = 0.524 for 25C.
#'@return the log10K at an ionic strength of zero
#' @export
log10K.IS.calculator = function(log10K.ref, I.ref, I, M.charge, X.charge, XM.charge, A = 0.524){
  y.XM = 10^(0.1*(XM.charge^2)*I.ref - (A*(XM.charge^2)*sqrt(I.ref))/(1+sqrt(I.ref)))
  y.M = 10^(0.1*(M.charge^2)*I.ref - (A*(M.charge^2)*sqrt(I.ref))/(1+sqrt(I.ref)))
  y.X = 10^(0.1*(X.charge^2)*I.ref - (A*(X.charge^2)*sqrt(I.ref))/(1+sqrt(I.ref)))
  log10K0 = log10K.ref - log10(y.XM/(y.M*y.X))
  y.XM = 10^(0.1*(XM.charge^2)*I - (A*(XM.charge^2)*sqrt(I))/(1+sqrt(I)))
  y.M = 10^(0.1*(M.charge^2)*I - (A*(M.charge^2)*sqrt(I))/(1+sqrt(I)))
  y.X = 10^(0.1*(X.charge^2)*I - (A*(X.charge^2)*sqrt(I))/(1+sqrt(I)))
  log10K = log10K.0 + log10(y.XM/(y.M*y.X))
}

#'Function that calculates the pKa of any Ionic strength using the Debye-Huckel equation
#'
#'@param pKa.ref pKa for the reference that you are referencing
#'@param I.ref Ionic strength of the reference constant
#'@param I ionic strength for the condition you want to calculate
#'@param X.charge charge of the base
#'@param HX.charge Charge of acid
#'@param A Debye-Huckel constant for the given temperature. Default A = 0.524 for 25C.
#'@return the log10K at an ionic strength of zero
#' @export
pKa.IS.calculator = function(pKa.ref, I.ref, I, X.charge, HX.charge, A = 0.524){
  y.HX = 10^(0.1*(HX.charge^2)*I.ref - (A*(HX.charge^2)*sqrt(I.ref))/(1+sqrt(I.ref)))
  y.H = 10^(0.1*I.ref - (A*sqrt(I.ref))/(1+sqrt(I.ref)))
  y.X = 10^(0.1*(X.charge^2)*I.ref - (A*(X.charge^2)*sqrt(I.ref))/(1+sqrt(I.ref)))
  pKa.0 = pKa.ref - log10(y.HX/(y.H*y.X))
  print("pKa.0")
  print(pKa.0)
  y.HX = 10^(0.1*(HX.charge^2)*I - (A*(HX.charge^2)*sqrt(I))/(1+sqrt(I)))
  y.H = 10^(0.1*I - (A*sqrt(I))/(1+sqrt(I)))
  y.X = 10^(0.1*(X.charge^2)*I - (A*(X.charge^2)*sqrt(I))/(1+sqrt(I)))
  pKa = pKa.0 + log10(y.HX/(y.H*y.X))
  print("pKa")
  print(pKa)
}
