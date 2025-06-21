# Soft-thresholding operator:

stoper=function(v,a){
  zz=((v-a)>0)*abs(v-a)-((-v-a)>0)*abs(-v-a)
  return(zz)
}





