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
