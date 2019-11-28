function angriff(p::BigInt,P::NfAbsOrdIdl,A::NfAbsOrdIdl,Mneu::NfAbsOrdIdl,B::NfAbsOrdIdl)
	k,a=quadratic_field(-p)
	c,mc=class_group(order(A))
	n=BigInt(order(c))
	println("n=",n)
	l=BigInt(order(preimage(mc,A)))
	println("l=",l)
	a=fmpz(n/l)
	println("a=",a)
	if Mneu==0
		return Hecke.power_class(P,a)==A
	else
		return Mneuneu=Hecke.reduce_ideal(Mneu*(Hecke.power_class(B,a))),l
	end
end	

function testmyangriff()
	m=BigInt(17)
	println("m=",m)
	p=BigInt(nextprime(1000))
	k,a=quadratic_field(-p)
	c,mc=class_group(k)
	println("c=",c)
	P=mc(c[1])
	println("p=",p)
	println("P=",P)
	a,A=diffiehellmanA(P)
	println("a=",a)
	println("A=",A)
	M=encode(p,m)
	println("M=",M)
	C,B=elgamalB(M,A,P)
	println("C=",C)
	Mneu=elgamalA(C,B,a)
	println("Mneu=",Mneu)
	Mneuneu,l=angriff(p,P,A,Mneu,B)
	println("Mneuneu=",Mneuneu)
	return M==Mneuneu
end
