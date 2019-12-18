#####
@doc Markdown.doc"""
attack(p::BigInt,P::NfAbsOrdIdl,A::NfAbsOrdIdl) -> BigInt

determines a s.th. P^a=A in O over Q[sqrt(-p)]
"""

function attack(p::BigInt,P::NfAbsOrdIdl,A::NfAbsOrdIdl)
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
	d=gcd(logP,BigInt(order(c))) 
	print("Fall ")
	if d==1
		println(1)
		return logA*(gcdx(logP,BigInt(order(c)))[2])
	else
		println(2)
		aneu=BigInt(logA/d)*gcdx(BigInt(logP/d),BigInt(BigInt(order(c))/d))[2]
		count=0
		Aneu=Hecke.power_class(P,fmpz(aneu))
		while isone(Hecke.reduce_ideal(A*inv(Aneu)))==false && count<d
			aneu+=BigInt(BigInt(order(c))/d)
			Aneu*=Hecke.power_class(P,fmpz(BigInt(order(c))/d))
			count+=1
		end
		println("count=",count)
		return aneu
	end
end	

function testmyattack(i::BigInt)
	p=BigInt(nextprime(i))
	println("p=",p)
	k,a=quadratic_field(-p)
	c,mc=class_group(k)
	println("c=",c)
	P=findnonprincipal(p,p*10)
	println("P=",P)
	a,A=diffiehellmanA(P)
	println("a=",a)
	println("A=",A)
	aneu=attack(p,P,A)
	println("aneu=",aneu)
	Aneu=Hecke.power_class(P,fmpz(aneu))
	return isone(Hecke.reduce_ideal(A*inv(Aneu)))
end
#####
@doc Markdown.doc"""
attack2(p::BigInt,P::NfAbsOrdIdl,A::NfAbsOrdIdl) -> BigInt

determines a s.th. P^a=A in O over Q[sqrt(-p)]
"""

function attack2(p::BigInt,P::NfAbsOrdIdl,A::NfAbsOrdIdl)
	k,a=quadratic_field(-p)
	c,mc=class_group(k)
	println("c=",c)
	logP=preimage(mc,P)
	println("logP=",logP)
	ordP=order(preimage(mc,P))
	println("ordP=",ordP)
	logA=preimage(mc,A)
	println("logA=",logA)
	h=hom(AbelianGroup([ordP]),c,[logP])
	X=haspreimage(h,logA)
	if X[1]
		return X[2][1]
	else
		print("fehler")
	end
end

function testmyattack2(i::BigInt)
	p=BigInt(nextprime(i))
	P=findnonprincipal(p,p*10)
	println("p=",p)
	println("P=",P)
	a,A=diffiehellmanA(P)
	println("a=",a)
	println("A=",A)
	aneu=attack2(p,P,A)
	println("aneu=",aneu)
	Aneu=Hecke.power_class(P,aneu)
	return isone(Hecke.reduce_ideal(A*inv(Aneu)))
end
#####
@doc Markdown.doc"""
bsgs(p::BigInt,P::NfAbsOrdIdl,A::NfAbsOrdIdl,C::NfAbsOrdIdl,B::NfAbsOrdIdl) -> BigInt

determines a s.th. P^a=A in O over Q[sqrt(-p)] via Babystep-Giantstep
"""

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
	#println("fertig1")
	for i=1:m
		push!(h,Hecke.reduce_ideal(h[i]*x))
		#println(i)
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
	aneu=bsgs(p,P,A,0*zk,0*zk)
	println("aneu=",aneu)
	return (a%order(c))==aneu
end
				
