function nucomp2(f1::QuadForm,f2:.QuadForm)
	D=abs(Hecke.discriminant(f1))
	L=root(divexact(D,4),4)
	if f1.a < f2.a
		x=f1
		f1=f2
		f2=x
	end
	u1=f1.a
	v1=f1.b
	w1=f1.c
	u2=f2.a
	v2=f2.b
	w2=f2.c
	s=divexact((v1+v2),2)
	m=v2-s
	(b,c,F)=gcdx(u2,u1)
	if mod(s,F)==0
		G=F
		Ax=G
		Bx=m*b
		By=divexact(u1,G)
		Cy=divexact(u2,G)
		Dy=divexact(s,G)
		bx=mod(Bx,By)
		by=By
		x=1
		y=0
		z=0
		while abs(by)>L && bx!=0
			q,t = divrem(by, bx)
			by=bx
			bx=t
			t=y-q*x
			y=x
			x=t
			z+=1	
		end
		if mod(z,2)==0
			by=-by
			y=-y
		end
		ax=G*x
		ay=G*y
		if z==0
			Q1=cx*bx
			cx=divexact(Q1-m,By)
			dx=
