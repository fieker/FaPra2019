function positive_definite_form(f::QuadForm)
	if f.a<0
		f.a=-f.a
		f.c=-f.c
	end
	return f
end

function isequal(f1::QuadForm,f2::QuadForm)
	f1=formreduce(f1)
	f2=formreduce(f2)
	return f1.a==f2.a && f1.b==f2.b && f1.c==f2.c
end

function formreduce(f::QuadForm)
	f=positive_definite_form(f)
	a = f.a
    b = f.b
    c = f.c
	if b<=(-a) || a<b
		q, r = divrem(b, 2*a)
		if r<0#??????????
			q=q-1
			r=r+2*a
		end
		if r>a
			r=r-2*a
			q=q+1
		end
		r=fmpz(r)
		q=fmpz(q)
		c=c-divexact((b+r)*q, 2)
		b=r
	end
	while a>c
		b=-b 
		x=a
		a=c
		c=x
        q, r = divrem(BigInt(b), BigInt(2*a))
		if r<0
			q=q-1
			r=r+2*a
		end
		r=fmpz(r)
		q=fmpz(q)
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
	return positive_definite_form(QuadForm(a, b, c))
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
	r=mod((y1*y2*n-x2*c2),v1)
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
	if mod(z,2)==1
		v2=-v2
		v3=-v3
	end
	return (z,d,v2,v3)
end

function nudupl(f::QuadForm)
	f=formreduce(f)
	a=f.a
	b=f.b
	c=f.c
	D=abs(Hecke.discriminant(f))
	L=root(divexact(D,4),4)#??
	(d1,u,v)=gcdx(b,a)
	A=divexact(a,d1)
	B=divexact(b,d1)
	C=mod((-c*u),A)
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
	if isequal(f1,f2)
		nudupl(f1)
	else
		f1=formreduce(f1)
		f2=formreduce(f2)
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
		#println("a1=",a1)
		#println("a2=",a2)
		(d,u,v)=gcdx(a2,a1)#?????????????
		#println("s=",s)
		#println("n=",n)
		#println("d=",d)
		#println("u=",u)
		#println("v=",v)
		if d==1
			#println("1.Fall")
			A=-u*n
			#println("A=",A)
			d1=d
			
			A=mod(A,a1)
			A1=a1-A
			if A1<A
				A=-A1
			end
			#println("A=",A)
			#println("d1=",d1)
			#println("A1=",A1)
			(z,d,v2,v3)=parteucl(a1,A,L)
			#println("z=",z)
			#println("d=",d)
			#println("v2=",v2)
			#println("v3=",v3)
		elseif mod(s,d)==0
			#println("2.Fall")
			A=-u*n
			d1=d
			a1=divexact(a1,d1)
			a2=divexact(a2,d1)
			s=divexact(s,d1)
			A=mod(A,a1)
			A1=a1-A
			if A1<A
				A=-A1
			end
			(z,d,v2,v3)=parteucl(a1,A,L)
		else
			#println("3.Fall")
			(d1,u1,v1)=gcdx(s,d) 
			if d1>1
				a1=divexact(a1,d1)
				a2=divexact(a2,d1)
				s=divexact(s,d1)
				d=divexact(d,d1)
			end
			c1=mod(c1,d)
			c2=mod(c2,d)
			l=mod((-u1*(u*c1+v*c2)),d)
			A=-u*divexact(n,d)+l*divexact(a1,d)
			A=mod(A,a1)
			A1=a1-A
			if A1<A
				A=-A1
			end
			(z,d,v2,v3)=parteucl(a1,A,L)
		end
		if z==0
			#println("special case")
			Q1=a2*v3
			Q2=Q1+n
			f=divexact(Q2,d)
			g=divexact((v3*s+c2),d)
			a3=d*a2
			c3=v3*f+g*d1#Fehler in Buch?????
			b3=2*Q1+b2
			f3=QuadForm(a3,b3,c3)
			#println("f3=",f3)	
			f1=i
			f2=h
			#println("Q1=",Q1)
			#println("Q2=",Q2)
			#println("f=",f)
			#println("g=",g)
			#println("a3=",a3)
			#println("c3=",c3)
			#println("b3=",b3)
			return formreduce(f3)
		end
		b=divexact((a2*d+n*v),a1)
		Q1=b*v3
		Q2=Q1+n
		f=divexact(Q2,d)
		e=divexact((s*d+c2*v),a1)
		Q3=e*v2
		Q4=Q3-s
		#println("Q4=",Q4)
		#println("v=",v)
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
end	

function formpowmod(f::QuadForm,n::fmpz)
	s=sign(n)
	n=abs(n)
	if n==0
		g=basistoform(fmpz[1 0;0 1],fundamental(Hecke.discriminant(f)))
	elseif n==1
		g=formreduce(f)
		#println(g)
	elseif n%2==0
		g=nudupl(formpowmod(f,divexact(n,fmpz(2))))
		#println(g)
	else
		g=nucomp(nudupl(formpowmod(f,divexact(n,fmpz(2)))),formreduce(f))
		#println(g)
	end
	if s==-1
		A=formtobasis(g)
		A[2]=-A[2]
		g=basistoform(A,fmpz(fundamental(Hecke.discriminant(g))))
	end
	return g
end

function simplecompose(f1::QuadForm,f2::QuadForm)
		D=Hecke.discriminant(f1)
		a1=f1.a
		a2=f2.a
		b1=f1.b
		b2=f2.b
		c1=f1.c
		c2=f2.c
		s=divexact(b1+b2,2)
		n=divexact(b1-b2,2)
		(d1,u1,v1)=gcdx(a1,a2)
		(d,u2,v2)=gcdx(d1,s)
		(u,v,w)=(u2*u1,u2*v1,v2)
		d0=gcd(d,gcd(c1,gcd(c2,n)))
		a3=divexact(d0*a1*a2,d^2)
		b3=b2+divexact(2*a2*(v*(s-b2)-w*c2),d)
		c3=divexact(b3^2-D,4*a3)
		f3=QuadForm(a3,b3,c3)
		return formreduce(f3)
end

function simplepower(f::QuadForm,n::fmpz)
	s=1
	if n <0
		s=-1
		n=fmpz(-1*n)
	end
	if n==0
		g=basistoform(fmpz[1 0;0 1],fundamental(Hecke.discriminant(f)))
	elseif n==1
		g=formreduce(f)
	elseif mod(n,2)==0
		g=simplecompose(simplepower(f,divexact(n,2)),simplepower(f,divexact(n,2)))
	else
		h=simplecompose(simplepower(f,divexact(n,2)),simplepower(f,divexact(n,2)))
		g=simplecompose(h,f)
	end
	if s==-1
		A=formtobasis(g)
		A[2]=-A[2]
		g=basistoform(A,fmpz(fundamental(Hecke.discriminant(g))))
	end
	return g
end
