import math

def fact(n):
    n = int(n)
    p = 1
    for i in range (2, n+1): 
        p *= i
    return p

# returns Clebsch-Gordan coefficient <a a1 b b1|c c1> = <j1 m1 j2 m2| JM>
# a, a1, b, b1, c, c1 = j1, m1, j2, m2, J, M
def CG(a, a1, b, b1, c, c1):
    
    if ( a1+b1!=c1 or a<abs(a1) or b<abs(b1) or c<abs(c1) 
        or a<0 or b<0 or c<0 or a+b<c or abs(a-b)>c or a-a1!=int(a-a1)
        or b-b1!=int(b-b1) ):
        return 0
    else:
        delta = math.sqrt(fact(a+b-c)*fact(a-b+c)*fact(-a+b+c)/fact(a+b+c+1))
        
        root = (math.sqrt(fact(a+a1)*fact(a-a1)*fact(c+c1)*fact(c-c1)*(2*c+1)
                 /fact(b+b1)/fact(b-b1)))
        z = 0
        s = 0
        while a+b-c1-z>=0 and b+c-a1-z>=0 and a-a1-z>=0 and c-c1-z>=0 and a+b+c+1-z>=0:
            
            s += (math.pow(-1,a-a1+z)*fact(a+b-c1-z)*fact(b+c-a1-z)
            /fact(z)/fact(a-a1-z)/fact(c-c1-z)/fact(a+b+c+1-z))
            z += 1
        clebsch = root*s/delta
        return clebsch

# returns 6j-symbol
# a, b, c, d, e, f = j1, j2, j3, j4, j5, j6 
def sixj(a, b, c, d, e, f):
   
    def triangle(a, b, c):
        if ((c<=a+b) and (abs(a-b)<=c)):
            return 1
        
    if ( (a,b,c,d,e,f>=0)
        and triangle(a,b,c)==1
        and triangle(b,d,f)==1
        and triangle(a,e,f)==1 
        and triangle(c,d,e)==1 
        and int(a+b-c)==a+b-c
        and int(b+d-f)==b+d-f
        and int(a+e-f)==a+e-f
        and int(c+d-e)==c+d-e ):
        
        def delta(a, b, c):
            return math.sqrt(fact(a+b-c)*fact(a-b+c)*fact(-a+b+c)/fact(a+b+c+1))
    
        p = math.pow(-1,a+c+d+f)*delta(a,b,c)*delta(b,d,f)/delta(a,e,f)/delta(c,d,e)
        n = 0
        s = 0
    
        while ( a-b+d+e-n>=0 and -b+c+e+f-n>=0 and a+c+d+f+1-n>=0 and a-b+c-n>=0 
                and -b+d+f-n>=0 and a+e+f+1-n>=0 and c+d+e+1-n>=0 ):
            
            s += ( math.pow(-1,n)*fact(a-b+d+e-n)*fact(-b+c+e+f-n)*fact(a+c+d+f+1-n)
            /fact(n)/fact(a-b+c-n)/fact(-b+d+f-n)/fact(a+e+f+1-n)/fact(c+d+e+1-n) )
            n += 1
        sixj = p*s
        return sixj
    else:
        return 0

# returns 9j-symbol
# a,b,c,d,e,f,g,h,j = j1,j2,j3,j4,j5,j6,j7,j8,j9
def ninej(a,b,c,d,e,f,g,h,j):
    
    xmin = max(abs(a-j), abs(b-f), abs(d-h))
    xmax = min(a+j, b+f, d+h) 
    s = 0
    x = xmin
    while x<=xmax:
        s += (math.pow(-1, 2*x)*(2*x+1)*sixj(a,b,c,f,j,x)
              *sixj(d,e,f,b,x,h)*sixj(g,h,j,x,a,d))
        x += 1
    return s
