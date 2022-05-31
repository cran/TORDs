
#' Function for generating the moment matrix and variance of the predicted response
#'
#' @param matrix Design matrix with the coefficients of the corresponding input factors
#'
#' @return The moment matrix and the prediction variance for a given design based on a third-order model
#'It gives unique prediction variance along with its frequencies.
#'@description This function generates the moment matrix and variance of the predicted response for a given design based on a third-order model, for measuring the rotatability of the design. The input should be the specified form of a design matrix with the coefficients of the corresponding input factors. A minimum number of centre points is to be used to ensure the non-singularity of X`X.
#' @export
#'@references
#'M. Hemavathi, Shashi Shekhar, Eldho Varghese, Seema Jaggi, Bikas Sinha & Nripes Kumar Mandal (2022)<DOI:10.1080/03610926.2021.1944213>." Theoretical developments in response surface designs: an informative review and further thoughts".
#' @examples
#' \dontrun{
#'library(TORDs)
#'library(TORDs)
#'Pred3.var(matrix)
#'}

Pred3.var<-function(matrix){
  v<-ncol(matrix)
  matfA<-matrix
  mat1<-matrix(,nrow=nrow(matfA),ncol=0)
  p=1
  while(p<=v){
    x1<-matrix(,nrow=nrow(matfA),ncol=0)
    x1<-(matfA[,p])^2
    mat1<-cbind(mat1,x1)
    p=p+1
  }
  ############normal matrix
  #matfA
  ##########sq terms mat1
  #mat1
  ##############interactions
  b1=1
  b2=1
  mat12<-matrix(,nrow=nrow(matfA),ncol=0)
  while(b1<v){
    b2=b1+1
    while(b2<=v){
      mat2<-matrix(,nrow=0,ncol=1)
      mat2<-matfA[,b1]*matfA[,b2]
      mat12<-cbind(mat12,mat2)
      b2=b2+1
    }

    b1=b1+1

  }
  #########interactions

  ############tricky portion
  new<-matfA

  o=1
  while(o<=v){
    matt<-matrix(,nrow=nrow(mat1),ncol=0)
    matt<-mat1[,c(o)]
    mat11<-mat1[,-c(o)]
    mat11<-cbind(matt,mat11)
    changed<-cbind(new[,o],mat11)
    mat123<-matrix(,nrow=nrow(new),ncol=0)
    q=1
    while(q<ncol(changed)){

      x1<-matrix(,nrow=nrow(new),ncol=0)
      x1<-changed[,1]*changed[,q+1]
      mat123<-cbind(mat123,x1)
      q=q+1
    }
    mat12<-cbind(mat12,new[,o],mat123)
    o=o+1
  }

  #################################checking


  #############3 factor interactions(wrong)
  l1=1
  l2=1
  l3=1
  int<-matrix(,nrow=nrow(mat1),ncol=0)
  while(l1<=v-2){
    l2=l1+1

    while(l2<=v-1){
      l3=l2+1
      while(l3<=v){
        mat2<-matrix(,nrow=0,ncol=1)
        mat2<-matfA[,l1]*matfA[,l2]*matfA[,l3]
        int<-cbind(int,mat2)
        l3=l3+1
      }
      l2=l2+1
    }

    l1=l1+1

  }
  #########################
  x_mat<-cbind(mat1,mat12,int)
  x_mat<-cbind((matrix(1,nrow=nrow(x_mat),ncol=1)),x_mat)
  rownames(x_mat)<-NULL
  colnames(x_mat)<-NULL
  x_prime_x<-t(x_mat)%*%x_mat


message("Moment matrix")
  N<-nrow(x_mat)
  mm<-x_prime_x/N
  print(mm)

  ###############
  x_matrix<-x_mat
  k1=1
  var<-c()
  while(k1<=nrow(x_matrix))
  {
    V=t(x_matrix[k1,])
    b<-t(V)
    v_y_hat<-V %*%solve(x_prime_x) %*% b
    var<-c(var,v_y_hat)
    k1<-k1+1
  }
  variance_of_esitmated_response<-round(var,digits = 4 )

  v2<-variance_of_esitmated_response
  v1<-unique(variance_of_esitmated_response)
  message("Predicted variance along with respective frequencies")
  for(j in v1){
    d<-length(which(v2==j))
    print(j)
    print(d)

  }

}
#}#################


