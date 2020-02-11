function diffiehellman0(m::fmpz)
	d=rand(m*fmpz(10)^3:m*fmpz(10)^6)
	while Hecke.isprime(d)==false || mod(-d,4)!=3
		d=rand(1000:m*fmpz(10)^6)
	end
	return -d
end

function find(d::fmpz)
	i=rand(1000:abs(d)*10)
	while Hecke.isprime(i)==false || jacobi(mod(d,i),i)!=1
		i=rand(1:fmpz(abs(d))*10)
	end
	I=[i 0;-lift(root(GF(i)(d),2)) 1]
	#print(I)
	return basistoform(I,d)
end

function encode(m::fmpz,d::fmpz)
	for i=m*fmpz(10)^8:(m*fmpz(10)^9-1)
		if Hecke.isprime(i) && jacobi(mod(d,i),i)==1
			M=[i 0;-lift(root(GF(i)(d),2)) 1]
			return basistoform(M,d)
		end
	end
end

function decode(M::QuadForm)
	M=formtobasis(M)[1]
	m=mod(M,10^9)
	return m
end

function diffiehellmanA(I::QuadForm,d::fmpz)
	a=rand(1000:fmpz(10)^6)
	A=formtobasis(formpowmod(I,fmpz(a)))
	#A=formtobasis(simplepower(I,fmpz(a)))
	while isone(A)############
		a=rand(1000:fmpz(10)^6)
		A=formtobasis(formpowmod(I,fmpz(a)))
		#A=formtobasis(simplepower(I,fmpz(a)))
	end
	#println("A=",A)
	return a,basistoform(A,d)
end

function diffiehellmanB(I::QuadForm,d::fmpz)
	b=rand(1000:fmpz(10)^6)
	B=formtobasis(formpowmod(I,fmpz(b)))
	while isone(B)
		b=rand(1000:fmpz(10)^6)
		B=formtobasis(formpowmod(I,fmpz(b)))
		#B=formtobasis(simplepower(I,fmpz(b)))
	end
	return b,basistoform(B,d)
end

function elgamalB(b::fmpz,A::QuadForm,M::QuadForm)
	K=formpowmod(A,b)
	#K=simplepower(A,b)
	C=nucomp(K,M)
	#C=simplecompose(K,M)
	return C
end

function elgamalA(a::fmpz,B::QuadForm,C::QuadForm)
	L=formpowmod(B,a)
	#L=simplepower(B,a)
	M=nucomp(L,C)
	#M=simplecompose(L,C)
	return M
end

function testmyelgamalwithforms(m::fmpz)
	d=diffiehellman0(m)
	println("d=",d)
	I=find(d)
	#println("I=",I)
	println(1)
	M=encode(m,d)
	#println("M=",M)
	println(2)
	a,A=diffiehellmanA(I,d)
	#println("a=",a)
	#println("A=",A)
	println(3)
	b,B=diffiehellmanB(I,d)
	#println("b=",b)
	#println("B=",B)
	println(4)
	C=elgamalB(b,A,M)
	#println("C=",C)
	println(5)
	Mneu=elgamalA(a,B,C)
	#println("Mneu=",Mneu)
	println(6)
	mneu=decode(Mneu)
	println("mneu=",mneu)
	println(7)
	return isequal(M,Mneu)
end

function smalltest()
	d=fmpz(-17)
	println("d=",d)
	m=fmpz(2)
	println("m=",m)
	I=fmpz[3 0;-lift(root(GF(3)(d),2)) 1]
	I=basistoform(I,d)
	println("I=",I)
	I=positive_definite_form(formreduce(positive_definite_form(I)))
	println("I=",I)
	M=fmpz[107 0;-lift(root(GF(107)(d),2)) 1]
	M=basistoform(M,d)
	println("M=",M)
	M=positive_definite_form(formreduce(positive_definite_form(M)))
	println("M=",M)
	a=fmpz(7)
	A=formpowmod(I,a)
	#A=simplepower(I,a)
	println("A=",A)
	A=positive_definite_form(formreduce(positive_definite_form(A)))
	println("A=",A)
	b=fmpz(9)
	B=formpowmod(I,b)
	#B=simplepower(I,b)
	println("B=",B)
	B=positive_definite_form(formreduce(positive_definite_form(B)))
	println("B=",B)
	C=elgamalB(b,A,M)
	println("C=",C)
	C=positive_definite_form(formreduce(positive_definite_form(C)))
	println("C=",C)
	Mneu=elgamalA(a,B,C)
	println("Mneu=",Mneu)
	Mneu=positive_definite_form(formreduce(positive_definite_form(Mneu)))
	println("Mneu=",Mneu)
	mneu=mod(formtobasis(Mneu)[1],100)
	println("mneu=",mneu)
	return m==mneu
end

#julia> smalltest()#mit nucomp, nudupl
#d=-17
#m=2
#I=<3, 2, 6>
#I=<3, 2, 6>
#M=<13, 20, 9>
#M=<2, 2, 9>
#A=<-27, 36, -6>
#A=<3, 6, -51>
#B=<3, 0, 5>
#B=<3, 0, 5>
#C=<6, -2, -55>
#C=<55, 112, 51>
#Mneu=<165, 12, 112>
#Mneu=<112, -12, 165>
#true

