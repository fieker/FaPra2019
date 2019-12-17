######
@doc Markdown.doc"""
findnonprincipal(p::BigInt,s::BigInt) -> NfAbsOrdIdl

returns a nonprincipal ideal in O over Q[sqrt(-p)] with first generator >= s (but only slightly)
"""
function findnonprincipal(p::BigInt,s::BigInt)
	k,a=quadratic_field(-p)
	#println(k)
	if s%10^8!=0
		@time max_order1(p)
	end
	@time zk=maximal_order(k)
	#println("zk bestimmt")
	i=nextprime(s)
	n=s*1000+1000
	while i<n
		 x=prime_decomposition(zk,fmpz(i))
		for j=1:length(x)
			if length(x)>1
				if isone(Hecke.reduce_ideal(x[j][1]))==false
					return x[j][1]
				end
			end
		end
		i=nextprime(i+1)
	end
	print("alle Ideale sind Hauptideale")
end

@doc Markdown.doc"""
max_order1(p::BigInt) 
sets maximal order of an AnticNumberField knowing that p is squarefree
"""
function max_order1(p::BigInt)
	k,a=quadratic_field(-p)
	if -p%4==1
		d=1//2*a+1//2
	else
		d=a
	end
	zk=Order(k,[k(1),d];check=false,cached=false)
	zk.ismaximal=1
	Hecke._set_maximal_order_of_nf(k,zk)
end
