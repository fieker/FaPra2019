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
	n=BigInt(order(c))
	println("n=",n)
	l=BigInt(order(preimage(mc,A)))
	println("l=",l)
	a=fmpz(n/l)
	println("a=",a)
	if C==0
		return Hecke.power_class(P,a)==A
	else
		return Mneuneu=Hecke.reduce_ideal(C*(Hecke.power_class(B,a))),l
	end

end	

function testmyangriff()
	m=BigInt(17)
	println("m=",m)
	p=BigInt(nextprime(10000000000))
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
	g=[]
	h=[A]
	d=0
	l=0
	x=Hecke.power_class(P,fmpz(-m))
	for i=1:m
		push!(g,Hecke.power_class(P,fmpz(i-1)))
	end
	for i=1:m
		push!(h,h[i]*x)
		#println(i)
		#println(h)
		for j=1:m
			println(Hecke.reduce_ideal(h[i+1]*inv(g[j])))
			if Hecke.reduce_ideal(h[i+1]*inv(g[j])).gen_one==1
				d=j
				l=i+1
				println()
				break
			end
		end
	end
	a=d+l*m
	#if C==0
		return a
	#else
	#	return Mneuneu=Hecke.reduce_ideal(C*(Hecke.power_class(B,fmpz(a))))
	#end
end

function testbsgs()
	m=BigInt(17)
	println("m=",m)
	p=BigInt(nextprime(10000000000000))
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
	a=bsgs(p,P,A,C,B)
	println("a=",a)
	#return Hecke.reduce_ideal(M*inv(Mneuneu))
end
				
