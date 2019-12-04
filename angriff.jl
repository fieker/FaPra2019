function encode_klein(p::BigInt,m::BigInt)
#include find.jl
	M=findnonprincipal(BigInt(p),m*10^3,m*10^4)
	return Hecke.reduce_ideal(M)
end

function decode_klein(M::NfAbsOrdIdl)
	m=divrem(M.gen_one,10^3)[1]
	return m
end

function angriff(p::BigInt,P::NfAbsOrdIdl,A::NfAbsOrdIdl,C::NfAbsOrdIdl,B::NfAbsOrdIdl)
	k,a=quadratic_field(-p)
	c,mc=class_group(order(A))
	n=collect(Primes.factor(BigInt(order(c))))
	ar =[]
	preP=BigInt(preimage(mc,P).coeff[1])
	preA=BigInt(preimage(mc,A).coeff[1])
	partP=[]
	partA=[]
	ar2=[]
	vgl=[]
	ar3=[]
	for i=1:length(n)
		push!(ar,n[i][1]^n[i][2])
		push!(partP,preP%ar[i])
		push!(partA,preA%ar[i])
		push!(ar2,div(ar[i],gcd(partP[i],ar[i])))
		push!(vgl,gcd(partP[i],ar[i])!=1)
		push!(ar3,partA[i]*gdcx(partP[i],ar2[i])[2])
	end
	crt(ar3,ar2)
end	

function testmyangriff()
	m=BigInt(17)
	println("m=",m)
	p=BigInt(nextprime(10))
	k,a=quadratic_field(-p)
	c,mc=class_group(k)
	println("c=",c)
	P=mc(c[1])
	println("p=",p)
	println("P=",P)
	a,A=diffiehellmanA(P)
	println("a=",a)
	println("A=",A)
	M=encode_klein(p,m)
	println("M=",M)
	C,B=elgamalB(M,A,P)
	println("C=",C)
	#Mneu=elgamalA(C,B,a)
	#println("Mneu=",Mneu)
	Mneuneu,l=angriff(p,P,A,C,B)
	println("Mneuneu=",Mneuneu)
	return Hecke.reduce_ideal(M*inv(Mneuneu))
end

function bsgs(p::BigInt,P::NfAbsOrdIdl,A::NfAbsOrdIdl,C::NfAbsOrdIdl,B::NfAbsOrdIdl)
	k,a=quadratic_field(-p)
	c,mc=class_group(order(A))
	n=BigInt(order(c))
	m=BigInt(ceil(sqrt(n)))
	println(m)
	g=[]
	h=[A]
	d=0
	l=0
	x=Hecke.power_class(P,fmpz(-m))
	for i=1:m
		push!(g,Hecke.power_class(P,fmpz(i-1)))
	end
	println("fertig1")
	for i=1:m
		push!(h,Hecke.reduce_ideal(h[i]*x))
		println(i)
		for j=1:m
			if isone(Hecke.reduce_ideal(h[i]*inv(g[j])))
				d=j-1
				l=i-1
				return d+l*m
			end
		end
	end
end

function testbsgs(p::Int64)
	#m=BigInt(17)
	#println("m=",m)
	p=BigInt(nextprime(p))
	k,a=quadratic_field(-p)
	zk=maximal_order(k)
	c,mc=class_group(k)
	println("c=",c)
	P=mc(c[1])
	println("p=",p)
	println("P=",P)
	a,A=diffiehellmanA(P)
	println("a=",a)
	println("A=",A)
	#M=encode_klein(p,m)
	#println("M=",M)
	#C,B=elgamalB(M,A,P)
	#println("C=",C)
	#Mneu=elgamalA(C,B,a)
	#println("Mneu=",Mneu)
	aneu=bsgs(p,P,A,0*zk,0*zk)
	println("aneu=",aneu)
	return (a%order(c))==aneu
	#return Hecke.reduce_ideal(M*inv(Mneuneu))
end
				
