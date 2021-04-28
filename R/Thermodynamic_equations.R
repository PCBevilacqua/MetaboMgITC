Wiseman.isotherm = function(V, X, M, H = 2534, K = 2880){
  x = X/M
  r = 1/(K*M)
  y = 0.5 + ((1 - x - r)/(2*sqrt(((1 + x + r)^2)-(4*x))))
  dQ.dX = H*V*y
}

Gibbs.free.energy = function(K, Temperature, R = 0.0019872){
  output = -R*(Temperature + 273.15)*log(K)
}

Entropy = function(dG, dH, Temperature){
  output = (dH - dG)/Temperature
}

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
