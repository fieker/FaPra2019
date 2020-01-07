function reduce(f::Array{Int64,1})
	a=f[1]
	b=f[2]
	c=f[3]
	#D=b*b-4*a*c
	while b<=(-a) || a<b
		r=Int64(b%(2+a))
		q=Int64((b-r)/(2+a))
		r=r-2*a
		q=q+1
		c=c-((b+r)*q)/2
		b=r
	end
	if a>c
		b=-b 
		x=a
		a=c
		c=x
		r=Int64(b%(2+a))
		q=Int64((b-r)/(2+a))
		r=r-2*a
		q=q+1
		c=c-((b+r)*q)/2
		b=r
	elseif a==c && b<0
		b=-b
	end
	return [a;b;c]
end

function compose(f::Array{Int64,1},g::Array{Int64,1})
	a1=f[1]
	b1=f[2]
	c1=f[3]
	a2=g[1]
	b2=g[2]
	c2=g[3]
	if a1>a2
		x=f
		f=g
		g=x
	end
	s=Int64((b1+b2)/2)
	n=b2-s
	if a2%a1==0
		y1=Int64(0)
		d=a1
	else
		y1=Int64(gcdx(a2,a1)[2])
	end
	if s%d==0
		y2=Int64(-1)
		x2=Int64(0)
		d1=d
	else
		x2=Int64(gcdx(s,d)[2])
		y2=Int64(-gcdx(s,d)[3])
	end
	v1=a1/d1
	v2=a2/d1
	r=(y1*y2*n-x2*c2)%v1
	b3=b2+2*v2*r
	a3=v1*v2
	c3=(c2*d1+r*(b2+v2*r))/v1
	return reduce([a3;b3;c3])
end

function parteucl(a::Int64,b::Int64,L::Int64)
	v=0
	d=a
	v2=1
	v3=b
	z=0
	while abs(v3)>L
		t2=v-gcdx(v3,t3)[2]*v2
		v=v2#d=v3
		v2=t2
		v3=t3
		z+=1
	end
	if z%2==0
		v2=-v2
		v3=-v3
	end
	return (z,d)
end

function nudupl(f::Array{Int64,1})
	a=f[1]
	b=f[2]
	c=f[3]
	D=b*b-4*a*c
	L=Int64(floor(sqrt(sqrt(D/4))))
	A=a/gcd(b,a)
	B=b/gcd(b,a)
	C=(-c*gcdx(b,a)[2])%A
	C1=A-C
	if C1<C
		C=-C1
	end
	(z,d)=parteucl(A,C,L)
	if z==0
		g=(B+v3+c)/d
		a2=d*d
		c2=v3*v3
		b2=b+(d+v3)^2-a2-c2
		c2=c2+g*gcd(b,a)
		return reduce([a2;b2;c2])
	end
	e=(c*v+B*d)/A
	g=(e*v2-B)/v
	b2=e*v2+v*g
	if d1>1
		b2=d1*b2
		v=d1*v
		v2=a1*v2
	end
	a2=d^2
	c2=v3^2
	b2=b2+(d+v3)^2-a2-c2
	a2=a2+e*v
	c2=c2+g*v2
	return reduce([a2;b2;c2])
end

function nucomp(f::Array{Int64,1},g::Array{Int64,1})
	a1=f[1]
	b1=f[2]
	c1=f[3]
	a2=g[1]
	b2=g[2]
	c2=g[3]
	D=b*b-4*a*c
	L=Int64(floor(sqrt(sqrt(D/4))))
	if a1<a2
		x=f
		f=g
		g=x
	end
	s(b1+b2)/2
	n=b2-s
	(d,u,v)=gcdx(a1,a2)
	if d==1
		A=-gcdx(a1,a2)[2]*n
		d1=1
		A=A%a1
		A1=a1-A
		if A1<A
			A=-A1
		end
		(z,d)=parteucl(a1,A)
	elseif s%d==0
		A=-gcdx(a1,a2)[2]*n
		d1=gcd(a1,a2)
		a1=a1/d1
		a2=a2/d1
		s=s/d1
		A=A%a1
		A1=a1-A
		if A1<A
			A=-A1
		end
		(z,d)=parteucl(a1,A)
	else
		(d1,u1,v1)=gcdx(s,d) 
		if d1>1
			a1=a1/d1
			a2=a2/d1
			s=s/d1
			d=d/d1
		end
		c1=c1%d
		c2=c2%d
		l=(-u1*(u*c1+v*c2))%d
		A=-u*(n/d)+l*(a1/d)
		A=A%a1
		if A1<A
			A=-A1
		end
		(z,d)=parteucl(a1,A)
	end
	if z==0
		Q1=a2*v3
		Q2=Q1+n
		h=Q2/d
		i=(v3*s+c2)/d
		a3=d*a2
		c3=v3*d+g*d1
		b3=2*Q1+b2
		return reduce([a3;b3;c3])
	end
	b=(a2*d+n*v)/a1
	Q1=b*v3
	Q2=Q1+n
	f=Q2/d
	e=(s*d+c2v)/a1
	Q3=e*v2
	Q4=Q3-s
	g=Q4/v
	if d1<1
		v2=d1*v2
		v=d1*v
	end
	a3=d*b+e*v
	c3=v3*f+g*v2
	b3=Q1+Q2+d1*(Q3+Q4)
	return reduce([a3;b3;c3])
end	
	
