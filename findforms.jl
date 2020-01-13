function findform(p::BigInt)
	g=[2,2,2]
	while gcd(g)!=1
		if rem(p,4)==1
			d=p
		else
			d=4*p
		end		
		b=fmpz(rand(-10000:10000))
		x=fmpz(d+b^2)
		while rem(x,4)!=0
			b=fmpz(rand(-10000:10000))
			x=fmpz(d+b^2)
		end
		d=fmpz(d)
		println("d=",d)
		y=div(x,fmpz(4))
		println("y=",y)
		if y<0
			sig(y)=-1
		else
			sig(y)=1
		end
		a=fmpz(rand(1:sig(y)*y))
		while rem(y,a)!=0
			a=fmpz(rand(1:sig(y)*y))
		end
		c=fmpz(div(y,a))
		g=[a,b,c]
	end
	f=QuadForm(g[1],g[2],g[3])
	return f
end

function findform2(p::BigInt)
	g=zeros(3)
	g[1]=rand(1:10000)
	g[2]=rand(-10000:10000)
	g[3]=rand(-10000:10000)
	f=QuadForm(g[1],g[2],g[3])
	i=1
	while Hecke.discriminant(f)!=p || gcd(g)!=1	
		println(i)	
		g[1]=rand(1:10000)
		g[2]=rand(-10000:10000)
		g[3]=rand(-10000:10000)
		f=QuadForm(g[1],g[2],g[3])
		i+=1
	end
	return f
end
		
#function findnontrivialform(p::BigInt)
#	f=formreduce(findform(p::BigInt))
#	while 
#		f=formreduce(findform(p::BigInt))
#	end
#	return f
#end
