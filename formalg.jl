function formreduce(f::QuadForm)
	#D=b^2-4*a*c
        a = f.a
        b = f.b
        c = f.c
	if b<=(-a) || a<b
		q, r = divrem(b, 2*a)
		if r>a
			r=r-2*a
			q=q+1
		end
		c=c-divexact((b+r)*q, 2)
		b=r
	end
	while a>c
		b=-b 
		x=a
		a=c
		c=x
        q, r = divrem(b, 2*a)
		if r>a
			r=r-2*a
			q=q+1
		end
		c=c-divexact((b+r)*q, 2)
		b=r
	end
	if a==c && b<0
		b=-b
	end
	return QuadForm(a, b, c)
end

#|b|<=a<=c,c>0 if a=|b| or a=c then b>0
function testmyformreduce(f::QuadForm)
	if Hecke.discriminant(f)<0 && f.a>0
		g=formreduce(f)
		a=g.a
		b=g.b
		c=g.c
		x=abs(b)<=a && a<=c
		if a==abs(b) || a==c
			x= x && b>0
		end
		return x
	else
		println("falscher input")
	end
end
#julia> J=findnonprincipal(BigInt(5),BigInt(50))
#<67, sqrt(-5)-14>
#Norm: 67
#Minimum: 67
#two normal wrt: 67
#julia> g=idealtoform(J,BigInt(5))
#my first form: <67, 28, 3>
#julia> testmyformreduce(g)
#false


#wird wahrscheinlich nicht benötigt#############################################
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
	return formreduce([a3;b3;c3])
end
################################################################################

function parteucl(a::fmpz,b::fmpz,L::fmpz)
	v=0
	d=a
	v2=1
	v3=b
	z=0
	while abs(v3)>L
		q, t3 =divrem(d,v3)
		t2=v-q*v2
		v=v2
		d=v3
		v2=t2
		v3=t3
		z+=1
	end
	if z%2==0
		v2=-v2
		v3=-v3
	end
	return (z,d,v2,v3)
end

function nudupl(f::QuadForm)	
	a=f.a
	b=f.b
	c=f.c
	D=abs(Hecke.discriminant(f))
	L=root(divexact(D,4),4)
	(d1,u,v)=gcdx(b,a)
	A=divexact(a,d1)
	B=divexact(b,d1)
	C=(-c*u)%A
	C1=A-C
	if C1<C
		C=-C1
	end
	(z,d,v2,v3)=parteucl(A,C,L)
	if z==0
		g=divexact((B*v3+c),d)
		a2=d^2
		c2=v3^2
		b2=b+(d+v3)^2-a2-c2
		c2=c2+g*d1
		return formreduce(QuadForm(a2,b2,c2))
	end
	e=divexact((c*v+B*d),A)
	g=divexact((e*v2-B),v)
	b2=e*v2+v*g
	if d1>1
		b2=d1*b2
		v=d1*v
		v2=d1*v2
	end
	a2=d^2
	c2=v3^2
	b2=b2+(d+v3)^2-a2-c2
	a2=a2+e*v
	c2=c2+g*v2
	return formreduce(QuadForm(a2,b2,c2))
end

function nucomp(f1::QuadForm,f2::QuadForm)
	i=f1
	h=f2
	if f1.a<f2.a
		x=f1
		f1=f2
		f2=x
	end
	a1=f1.a
	b1=f1.b
	c1=f1.c
	a2=f2.a
	b2=f2.b
	c2=f2.c
	D=abs(Hecke.discriminant(f1))
	L=root(divexact(D,4),4)
	s=divexact(b1+b2,2)
	n=b2-s
	(d,u,v)=gcdx(a1,a2)#u=0 | v=0 ??
	if d==1
		A=-u*n
		d1=d
		A=A%a1
		A1=a1-A
		if A1<A
			A=-A1
		end
		println("a1=",a1)
		println("A=",A)
		println("L=",L)
		(z,d,v2,v3)=parteucl(a1,A,L)
	elseif s%d==0
		A=-u*n
		d1=d
		a1=divexact(a1,d1)
		a2=divexact(a2,d1)
		s=divexact(s,d1)
		A=A%a1
		A1=a1-A
		if A1<A
			A=-A1
		end
		(z,d,v2,v3)=parteucl(a1,A,L)
	else
		(d1,u1,v1)=gcdx(s,d) 
		if d1>1
			a1=divexact(a1,d1)
			a2=divexact(a2,d1)
			s=divexact(s,d1)
			d=divexact(d,d1)
		end
		c1=c1%d
		c2=c2%d
		l=(-u1*(u*c1+v*c2))%d
		A=-u*divexact(n,d)+l*divexact(a1,d)
		A=A%a1
		A1=a1-A
		if A1<A
			A=-A1
		end
		(z,d,v2,v3)=parteucl(a1,A,L)
	end
	if z==0
		Q1=a2*v3
		Q2=Q1+n
		f=div(Q2,d)
		g=div((v3*s+c2),d)
		a3=d*a2
		c3=v3*d+g*d1
		b3=2*Q1+b2
		f3=QuadForm(a3,b3,c3)	
		f1=i
		f2=h
		return f3#######
	end
	b=divexact((a2*d+n*v),a1)
	Q1=b*v3
	Q2=Q1+n
	f=divexact(Q2,d)
	e=divexact((s*d+c2*v),a1)
	Q3=e*v2
	Q4=Q3-s
	println("Q4=",Q4)
	println("v=",v)
	g=divexact(Q4,v)
	if d1>1
		v2=d1*v2
		v=d1*v
	end
	a3=d*b+e*v
	c3=v3*f+g*v2
	b3=Q1+Q2+d1*(Q3+Q4)
	f3=QuadForm(a3,b3,c3)
	f1=i
	f2=h
	return formreduce(f3)
end	
	
function formpowmod(f::QuadForm,n::fmpz)
	if n==1
		g=f
	elseif n%2==0
		g=nudupl(formpowmod(f,divexact(n,fmpz(2))))
	else
		g=nucomp(nudupl(formpowmod(f,divexact(n,fmpz(2)))),f)
	end
	return g
end
