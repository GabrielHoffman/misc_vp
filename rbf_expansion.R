

# RBF features map basis expansion 
# Jan 8, 2025
# Kernel method are quadratic in the number of samples
# since the Gram matrix is n x n
# Instead, create the explicit feature mapping using rbf_basis()
# and apply it to all data points with rbf_basis_expand()
# Now methods are linear in n and quadratic in the basis size
# But currently, rbf_basis() in only reproduces rdf() for scalar inputs
# How to expand this to multidimension observations?
# See https://stats.stackexchange.com/questions/69759/feature-map-for-the-gaussian-kernel
# and Eqn 9 here: https://arxiv.org/pdf/1109.4603

# For extension to kNN see https://davpinto.github.io/fastknn/

# Create kernal value for two points
rbf = function(x1, x2, sigmaSq){

	# sigmaSq * exp(- sum((x1-x2)^2) / (2*l^2))

	# https://andrewcharlesjones.github.io/journal/rbf.html
	exp(- sum((x1-x2)^2) / (2*sigmaSq))
}

# basis expansion for 1 point
rbf_basis = function(x, sigmaSq, K){

	# res = lapply(seq(K), function(k){
		# exp(-x^2/(2*l^2)) * x^k / (l^k * sqrt(factorial(k)))
		# exp(-1/(2*l^2)*x^2) / l^k / sqrt(factorial(k)) * x^k
		# })
	# sqrt(sigmaSq) * do.call(rbind, res)

	# https://stats.stackexchange.com/questions/69759/feature-map-for-the-gaussian-kernel
	# https://www.csie.ntu.edu.tw/~cjlin/talks/kuleuven_svm.pdf
	A = exp(-x^2/(2*sigmaSq))

	if( K > 0 ){
		res = sapply(seq(K), function(k){
			sqrt(1/(factorial(k)*sigmaSq^k)) * x^k
			})

		res = matrix(c(A, A*res), ncol=1)
	}else{
		res = matrix(A, ncol=1)
	}
	res
}

# Rcpp implementation
sourceCpp(code = "
#include <cmath>
#include<Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector rbf_basis1(const double &x, const double &sigmaSq, const int &K){
	double a = exp(-1.0*pow(x, 2) / (2.0*sigmaSq));

	NumericVector b(K+1);
	b[0] = 1;
	for(int k=1; k<=K; k++){
		b[k] = sqrt(1/(tgamma(k+1) * pow(sigmaSq,k))) * pow(x, k);
	}

	return a * b;
}

// [[Rcpp::export]]
NumericVector rbf_basis_expand(const NumericVector &x, const double &sigmaSq, const int &K){

	NumericMatrix res(x.size(), K+1);

	// set values for each data point
	for(int i=0; i<x.size(); i++){
		res(i,_) = rbf_basis1(x[i], sigmaSq, K);
	}

	// name columns
	CharacterVector cn(K+1);
	for(int k=0; k<=K; k++){
		cn(k) = \"b\" + std::to_string(k);
	}
	colnames(res) = cn;

	return res;
}
")

rbf_basis_expand(seq(4), 1, K)



.1
K = 3
rbf_basis(x, 1, K)
rbf_basis1(x, 1, K)


rbf_basis_expand(seq(4), 1, K)



n = 2
x = rnorm(n)
y = rnorm(n)

rbf(x,y, 1)
rbf(y,x, 1)

X1 = rbf_basis(x, 1, 8)
X2 = rbf_basis(y, 1, 8)

crossprod(X1, X2)



sigSq = 0.1
N = 100
df = expand.grid(i = seq(N), 
				j = seq(N))
df$x = (df$i - N/2)/(N/2)
df$y = (df$j - N/2)/(N/2)

with(df, plot(x,y, col=(x-y)^2))

# not a valid kernel
ggplot(df, aes(x,y, fill=(x-y)^2)) +
	geom_tile() +
	theme_classic() + 
	coord_fixed(expand=FALSE) +
	scale_fill_gradient2(low="white", high="red") + 
	ggtitle("Euclidian")

sigmaSq = .1
ggplot(df, aes(x,y, fill=exp(-1 / (2*sigmaSq) * (x-y)^2))) +
	geom_tile() +
	theme_classic() + 
	coord_fixed(expand=FALSE) +
	scale_fill_gradient2(low="white", high="red") + 
	ggtitle("Gaussian")




C = matrix(NA, N, N)
for(idx in seq(nrow(df))){
	C[df$i[idx], df$j[idx]] = rbf(df$x[idx], df$y[idx], sigSq )
}
image(C)

K = 4
G = matrix(NA, N, N)
for(idx in seq(nrow(df))){

	x1 = rbf_basis(df$x[idx], sigSq, K)
	x2 = rbf_basis(df$y[idx], sigSq, K)

	G[df$i[idx], df$j[idx]] = crossprod(x1, x2)
}
image(G)

K = 14
X = lapply(unique(df$x), function(x) rbf_basis(x, sigSq, K) )
X = do.call(cbind, X)

image(crossprod(X))

plot(C, crossprod(X))
abline(0, 1, col="red")

K = 89
X = rbf_basis_expand(unique(df$x), sigSq, K)
image(tcrossprod(X))

# Saturation curve of K depends strongly on sigSq
sigSq = 1e-2
k_values = seq(0, 100)
v = sapply(k_values, function(K){
	X = rbf_basis_expand(unique(df$x), sigSq, K)
	sum(colVars(X))
	})
plot(k_values, v)

# Test with kNN
#
# Compute kNN, and compute kernel matrix
library(cellpypes)
library(Matrix)
library(irlba)
n = 200
p = 30
# fmat <- matrix(rnorm(n*p), ncol=n)
fmat = t(seq(-10, 10, length.out=n))
nn <- find_knn(fmat,k=100)

C = matrix(NA, n, n)
for(i in seq(n)){
	C[i,nn$idx[i,]] = nn$dist[i,]
	C[nn$idx[i,],i] = nn$dist[i,]
}
C = max(C, na.rm=TRUE) - C
C = C / max(C, na.rm=TRUE)
C[is.na(C)] = 0
image(C)

plot(eigen(C)$values)


# Convert kNN similarities
# to features
C = Matrix(0, nrow=n, ncol=n)
for(i in seq(n)){
	C[i,nn$idx[i,]] = nn$dist[i,]
	C[nn$idx[i,],i] = nn$dist[i,]
}
m = max(C, na.rm=TRUE) 
C@x = m - C@x
C@x = C@x / max(C@x)
diag(C) = 1

image(as.matrix(C))

dcmp = partial_eigen(C, 20)

plot(dcmp$values)








