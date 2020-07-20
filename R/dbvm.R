foo<-function(y, L, prior){
  mvnfast::dmvn(0, 0, 1)
  CholWishart::rInvCholWishart(1, 4, diag(2))[, , 1]
}
