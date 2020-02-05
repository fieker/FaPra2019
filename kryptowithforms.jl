function diffiehellman0(m::fmpz)
	d=rand(1:m*10)
	while Hecke.isprime(d)==false || mod(-d,4)!=1
		d=rand(1:m*10)
	end
	return -d
end

function find(d::fmpz)
	i=rand(1:abs(d)*10)
	while Hecke.isprime(i)==false || jacobi(mod(d,i),i)!=1
		i=rand(1:abs(d)*10)
	end
	I=[i 0;lift(root(GF(i)(d),2)) 1]
	print(I)
	return basistoform(I,d)
end

function encrypt(m::fmpz,d::fmpz)
	for i=m*10:(m*10^2-1)
		if Hecke.isprime(i) && jacobi(mod(d,i),i)==1
			M=[i 0;lift(root(GF(i)(d),2)) 1]
			return basistoform(M,d)
		end
	end
end

function decrypt(M::QuadForm)
	M=formtobasis(M)[1]
	m=div(M,10)
	return m
end

function diffiehellmanA(I::QuadForm,d::fmpz)
	a=rand(10:fmpz(10)^2)
	A=formtobasis(formpowmod(I,fmpz(a)))
	while isone(A)
		a=rand(10:fmpz(10)^2)
		A=formtobasis(formpowmod(I,fmpz(a)))
	end
	return a,basistoform(A,d)
end

function diffiehellmanB(I::QuadForm,d::fmpz)
	b=rand(10:fmpz(10)^2)
	B=formtobasis(formpowmod(I,fmpz(b)))
	while isone(B)
		b=rand(10:fmpz(10)^2)
		B=formtobasis(formpowmod(I,fmpz(b)))
	end
	return b,basistoform(B,d)
end

function elgamalB(b::fmpz,A::QuadForm,M::QuadForm)
	K=formpowmod(A,b)
	C=nucomp(K,M)
	return C
end

function elgamalA(a::fmpz,B::QuadForm,C::QuadForm)
	L=formpowmod(B,a)
	M=nucomp(L,C)
	return M
end

function testmyelgamalwithforms(m::fmpz)
	d=diffiehellman0(m)
	println("d=",d)
	I=find(d)
	println("I=",I)
	M=encrypt(m,d)
	println("M=",M)
	a,A=diffiehellmanA(I,d)
	println("a=",a)
	println("A=",A)
	b,B=diffiehellmanB(I,d)
	println("b=",b)
	println("B=",B)
	C=elgamalB(b,A,M)
	println("C=",C)
	Mneu=elgamalA(a,B,C)
	println("Mneu=",Mneu)
	mneu=decrypt(Mneu)
	println("mneu=",mneu)
	return mneu==m
end
	

