function encode_klein(p::BigInt,m::BigInt)
#include find.jl
	M=findnonprincipal(BigInt(p),m*10^3,m*10^4)
	return Hecke.reduce_ideal(M)
end

function decode_klein(M::NfAbsOrdIdl)
	m=divrem(M.gen_one,10^3)[1]
	return m
end

function angriff(p::BigInt,P::NfAbsOrdIdl,A::NfAbsOrdIdl)
	#k,a=quadratic_field(-p)
	c,mc=class_group(order(A))
	n=collect(Primes.factor(BigInt(order(c))))
	#Primfaktorzerlegung der Ordnung der Klassengruppe
	ar =[]
	#teilerfremde Komponenten von n
	logP=BigInt(preimage(mc,P).coeff[1])
	println("logP=",logP)
	println("ordP=",order(preimage(mc,P)))
	logA=BigInt(preimage(mc,A).coeff[1])
	println("logA=",logA)
	partP=[]
	#Komponenten in den Zni
	partA=[]
	ar2=[]
	#partP soll teilerfremd zu den Einträgen in ar sein deshalb durch den ggt teilen
	vgl=[]
	#vergleiche wo ar2 und ar unterschiedlich sind
	ar3=[]
	#lag(A)/log(P) 
	if gcd(logP,BigInt(order(c)))==1
		println(1)
		return logA*(gcdx(logP,BigInt(order(c)))[2]),0
	end
	for i=1:length(n)
		push!(ar,n[i][1]^n[i][2])
		push!(partP,logP%ar[i])
		push!(partA,logA%ar[i])
		if fmpz(div(ar[i],gcd(partP[i],ar[i])))!=1
			push!(ar2,fmpz(div(ar[i],gcd(partP[i],ar[i]))))
		end
	end
	for i=1:length(ar2)
		push!(vgl,gcd(partP[i],ar[i])!=1)
		push!(ar3,fmpz(partA[i]*gcdx(div(partP[i],BigInt(gcd(partP[i],ar2[i]))),BigInt(ar2[i]))[2]))
	end
	println(ar)
	println(ar2)
	println(vgl)
	x=1
	#for i=1:length(vgl)
	#	if vgl[i]==true
	#		x*=(ar[i])/BigInt(ar2[i])
	#	end
	#end
	if length(ar2)>1
		println(2)
		println(logA*(gcdx(BigInt(logP/BigInt(gcd(logP,order(c)))),BigInt(order(c)))[2]))
		return crt(ar3,ar2),x
	else
		println(3)
		println(logA*(gcdx(BigInt(logP/BigInt(gcd(logP,order(c)))),BigInt(order(c)))[2]))
		return ar3[1],x
	end
end	

function testmyangriff(i::BigInt)
	p=BigInt(nextprime(i))
	k,a=quadratic_field(-p)
	c,mc=class_group(k)
	println("c=",c)
	#P=mc(c[1])
	P=findnonprincipal(p,p*10,p*200)
	println("p=",p)
	#println("P=",P)
	a,A=diffiehellmanA(P)
	println("a=",a)
	#println("A=",A)
	aneu,x=angriff(p,P,A)
	println("aneu=",aneu)
	println("x=",x)
	Aneu=Hecke.power_class(P,fmpz(aneu))
	count=1
	while isone(Hecke.reduce_ideal(A*inv(Aneu)))==false && count<order(c)
    	Aneu=Hecke.reduce_ideal(Aneu*Hecke.power_class(P,fmpz(x)))
		print(count,"  ")
		count+=1
    end
	return isone(Hecke.reduce_ideal(A*inv(Aneu)))
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
				
