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
			r=r+abs(2*a)
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
			q=q-sign(2*a)
			r=r+abs(2*a)
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

function parteucl(a::fmpz,b::fmpz,L::fmpz)
	v=0
	d=a
	v2=1
	v3=b
	z=0
	while abs(v3)>L	
		q, t3 =divrem(d,v3)
		if t3<0
			t3=t3+abs(v3)
			q=q-sign(v3)
		end
		t2=v-q*v2
		v=v2
		d=v3
		v2=t2
		v3=t3
		z+=1
	end
	if iseven(z)==false
		v2=-v2
		v3=-v3
	end
	return (z,d,v2,v3,v)
end

function nudupl(f::QuadForm)
	f=formreduce(f)
	a=f.a
	b=f.b
	c=f.c
	D=abs(Hecke.discriminant(f))
	L=root(divexact(D,4),4)
	(d1,u,v)=gcdx(b,a)#v wird nicht gebraucht,invmod falls gcd(b,a)=1?
	A=divexact(a,d1)
	B=divexact(b,d1)
	C=mod((-c*u),A)
	C1=A-C
	if C1<C
		C=-C1
	end
	(z,d,v2,v3,v)=parteucl(A,C,L)
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
	e=QuadForm(1,0,fundamental(Hecke.discriminant(f1)))
	if isequal(f1,e)
		return f2
	elseif isequal(f2,e)
		return f1
	elseif isequal(f1,f2)
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
		(d,u,v)=gcdx(a2,a1)#v wird nur in Fall3 gebraucht
		if d==1
			A=-u*n
			d1=d
			
			A=mod(A,a1)
			A1=a1-A
			if A1<A
				A=-A1
			end
			(z,d,v2,v3,v)=parteucl(a1,A,L)
		elseif mod(s,d)==0
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
			(z,d,v2,v3,v)=parteucl(a1,A,L)
		else
			(d1,u1,v1)=gcdx(s,d)#v1 wird nicht gebraucht
			if d1>1
				a1=divexact(a1,d1)
				a2=divexact(a2,d1)
				s=divexact(s,d1)
				d=divexact(d,d1)
			end
			c_1=mod(c1,d)
			c_2=mod(c2,d)
			l=mod((-u1*(u*c_1+v*c_2)),d)
			A=-u*divexact(n,d)+l*divexact(a1,d)
			A=mod(A,a1)
			A1=a1-A
			if A1<A
				A=-A1
			end
			(z,d,v2,v3,v)=parteucl(a1,A,L)
		end
		if z==0
			Q1=a2*v3
			Q2=Q1+n
			f=divexact(Q2,d)
			g=divexact((v3*s+c2),d)
			a3=d*a2
			c3=v3*f+g*d1#Fehler in Buch?????
			b3=2*Q1+b2
			f3=QuadForm(a3,b3,c3)
			f1=i
			f2=h
			return formreduce(f3)
		end
		b=divexact((a2*d+n*v),a1)
		Q1=b*v3
		Q2=Q1+n
		f=divexact(Q2,d)
		e=divexact((s*d+c2*v),a1)
		Q3=e*v2
		Q4=Q3-s
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
	elseif iseven(n)
		g=nudupl(formpowmod(f,divexact(n,fmpz(2))))
	else
		g=nucomp(nudupl(formpowmod(f,divexact(n,fmpz(2)))),formreduce(f))
	end
	if s==-1
		A=formtobasis(g)
		A[2]=-A[2]
		g=basistoform(A,fmpz(fundamental(Hecke.discriminant(g))))
	end
	return g
end

function formpowmodbit(f::QuadForm,n::fmpz)
	s=sign(n)
	n=abs(n)
	if n==0
		g=basistoform(fmpz[1 0;0 1],fundamental(Hecke.discriminant(f)))
	elseif n==1
		g=formreduce(f)
	else
		f=formreduce(f)
		g=basistoform(fmpz[1 0;0 1],fundamental(Hecke.discriminant(f)))
		j=1
		for i in bits(n)
			if i
				g=nucomp(g,f)
			end
			if j<length(bits(n))
				g=nudupl(g)
				j+=1
			end
		end
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
	elseif iseven(n)
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
