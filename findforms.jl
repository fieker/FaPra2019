function findform(p::BigInt)
	g=[2,2,2]
	while gcd(g)!=1
		if rem(p,4)==1
			d=p
		else
			d=4*p
		end
		b=fmpz(rand(-floor(sqrt(d)):floor(sqrt(d))))
		x=fmpz(d-b^2)
		while rem(x,4)!=0
			b=fmpz(rand(-floor(sqrt(d)):floor(sqrt(d))))
			x=fmpz(d-b^2)
		end
		d=fmpz(d)
		y=div(x,fmpz(4))
		a=fmpz(rand(1:y))
		while rem(y,a)!=0
			a=fmpz(rand(1:y))
		end
		c=fmpz(div(y,fmpz(4)))
		g=[a,b,c]
	end
	f=QuadForm(g[1],g[2],g[3])
	return f
end

#function findnontrivialform(p::BigInt)
#	f=formreduce(findform(p::BigInt))
#	while 
#		f=formreduce(findform(p::BigInt))
#	end
#	return f
#end
