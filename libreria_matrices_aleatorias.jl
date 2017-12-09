
__precompile__() #Se precompila el paquete.

module herramienta

export GOE

"""GOE(dim)

GOE es una función que realiza el ensamble de matrices aleatorias GOE cuya entrada corresponde a la dimensión."""

function GOE(dim) #Se crea la función para en ensamble de matrices aleatorias GOE
   (x -> (x+x')/2)(randn(dim,dim)) #Se toma la definición de matrices GOE
end

export GUE

"""GUE(dim)

GUE es una función que realiza el ensamble de matrices aleatorias GUE, cuya entrada corresponda a la dimensión."""

function GUE(dim) #Se crea la función para en ensamble de matrices aleatorias GUE
  (x -> (x+x')/2)(randn(dim,dim)+im*randn(dim,dim)) #Se toma de la definición de matrices GUE
end

export random_state

"""random_state(dim)

random_state es una función que crea estados aleatorios factorizables para ello se factoriza al estado con un 
espín (1,0) y un estado aleatorio. Esta función toma como entrada la dimensión y a la salida entrega un 
estado factorizado y normalizado."""

function random_state(dim)
    v=kron([1,0],randn(dim,1)) #El estado se puede factorizar con el spin (1,0), y un estado aleatorio
    v=v/norm(v) #Se normaliza el estado
    return v #Se regresa el estado normalizado
end

export projector

"""projector(state)

projector es una función que aplica el operador de proyección, su entrada es un estado cualquiera."""

function projector(state) #Se crea el operador de proyección
    return state*state'
end

export partial_trace_pure_bipartite_mat

"""partial_trace_pure_bipartite_mat(state,dim,system) 

partial_trace_pure_bipartite_mat es una función que calcula la traza parcial, sus entradas corresponden 
al estado, la dimensión y el sistema."""

#En esta función se calcula la traza parcial 
function partial_trace_pure_bipartite_mat(state,dim,system)
    
    dimtotal = length(state)[1] #Tomamos la dimensión de la matriz de densidad 
    
    dimcomp = Int(dimtotal/dim) #Se calcula la dimensión compuesta como la dimensión total entre la 
    #dimensión del sistema respecto al cual se tomará la traza parcial
    
    psi = reshape(state,(dimcomp,dim))'#Se define una nueva matriz con las dimensiones del sistema
    
    if system == 1#Si el sistema es igual al sistema 1
    
        psi=conj(psi) #Se calcula el conjugado
        
        return psi*psi' #Se regresa el proyector en ese estado
        
        elseif system == 2 #Si no, se debe calcular el producto interno
        
        return psi'*psi #Se regresa el produnto interno
    
    end

end

export H

"""H(a,n,g)

H es una función que calcula la pureza de estados generados aleatoriamente con una distribución 
gaussiana a distintos tiempos. Toma como entradas dos parámetros (a y g) y la dimensión (n) del hamiltoniano
del ambiente libre de interacciones"""

function H(a,n,g)
    
    Sz = [1 0;0 -1] #Se crea el operador matriz de pauli Sigma Z.
    
    goe = herramienta.GOE(2*n) #Definimos el hamiltoniano de interacción entre el sistema y el ambiente_1.
    
    gue = herramienta.GUE(n) #Definimos el hamiltoniano del ambiente libre de interacciones.
    
    h = (a/2)*kron(Sz,eye(n))+ kron(eye(2),gue)+ g*goe #Se define el hamiltoniano completo 
    #(incluido el sistema central).
    
    t = 0:0.1:100 #Se crea una lista de tiempo.
    
    V = herramienta.random_state(n) #Se crea un estado aleatorio de dimensión n.
    
    rho = [] #Se crea una lista para guardar las matrices de densidad.
    
    for i in 1:length(t) #En este for se evoluciona el estado dado en V usando el operador de evolución temporal
      
        psi = (expm(-im*h*t[i]))*(herramienta.projector(V))*(expm(im*h*t[i])) #Se calcula la matriz de densidad
        #al tiempo t del estado V.
       
        
        push!(rho,psi) #Se guarda en una lista la evolución de la matriz de densidad a distintos tiempos.
        
    end
    
    partial = [] #Creamos una lista para guardar la traza parcial.
    
    for i in 1:length(rho) #En este for se calcula la traza parcial para la matriz de densidad a cada tiempo
        
        trace = herramienta.partial_trace_pure_bipartite_mat(rho[i],2,1)
        
        push!(partial,trace) #Se guarda la traza parcial en una lista.
    
    end
    
    densit2 = [] #Se crea una lista para guardar la matriz al cuadrado.
    
    for i in 1:length(rho) #En este for se calcula el cuadrado de la traza parcial a cada tiempo t de la lista.
    
        densit = (partial[i])*(partial[i]) #Se eleva al cuadrado.
        
        push!(densit2,densit) #Se guarda la matriz obtenida.

    end
    
    traza = [] #Finalmente se calcula la traza de las matrices obtenidas en el paso anterior.
    
    for i in 1:length(densit2)
        
        traz = trace(densit2[i]) #Se calcula la traza para cada elemento de la lista.
        
        push!(traza,traz) #Guardamos el valor de la traza a cada tiempo t.
    
    end
    
    return t, traza, gue, goe #Se regresa la lista de tiempo, la lista con la traza al tiempo t, y 
    #las matrices GUE (hamiltoniano libre del ambiente) y la matriz GOE(hamiltoniano de interacción).
    
    
end

export H2

"""H2(a,n,g,G)

H2 es una función que calcula la pureza de estados generados aleatoriamente con una distribución 
gaussiana a distintos tiempos, considerando dos ambientes. Toma como entradas tres parámetros (a,g,G) 
y la dimensión (n) del hamiltoniano del ambiente libre de interacciones"""

function H2(a,n,g,G)
    
    Sz = [1 0;0 -1] #Se crea el operador matriz de pauli Sigma Z.
    
    goe = herramienta.GOE(2*n) #Definimos el hamiltoniano de interacción entre el sistema y el ambiente_1.
    
    gue = herramienta.GUE(n) #Definimos el hamiltoniano del ambiente libre de interacciones.
    
    t = 0:0.1:100 #Se crea una lista de tiempo.
    
    H = (a/2)*kron(Sz,eye(n),eye(n))+kron(eye(2),gue,eye(n))+g*kron(goe,eye(n))+kron(eye(2),eye(n),herramienta.GUE(n))+G*kron(eye(2),herramienta.GUE(n^2))
 
    V = herramienta.random_state(n*n) #Se crea un estado aleatorio de dimensión nxn.
    
    rho = [] #Se crea una lista para guardar las matrices de densidad.
    
    for i in 1:length(t) #En este for se evoluciona el estado dado en V usando el operador de evolución temporal.
    
        psi = (expm(-im*H*t[i]))*(herramienta.projector(V))*(expm(im*H*t[i])) #Se calcula la matriz de densidad 
        #al tiempo t del estado V.

   
        push!(rho,psi) #Se guarda en una lista la evolución de la matriz de densidad a distintos tiempos.
    
    end

    partial = [] #Creamos una lista para guardar la traza parcial.

    for i in 1:length(rho)

        trace = herramienta.partial_trace_pure_bipartite_mat(rho[i],2,1) #En este for se calcula la traza 
        #parcial para la matriz de densidad a cada tiempo.
        
        push!(partial,trace) #Se guarda la traza parcial en una lista.

    end
    
    densit2 = [] #Se crea una lista para guardar la matriz al cuadrado.
    
    for i in 1:length(rho) #En este for se calcula el cuadrado de la traza parcial a cada tiempo t de la lista.

        densit = (partial[i])*(partial[i]) #Se eleva al cuadrado.
        
        push!(densit2,densit) #Se guarda la matriz obtenida.

    end
    
    traza = [] #Finalmente se calcula la traza de las matrices obtenidas en el paso anterior.
    
    for i in 1:length(densit2) #Se calcula la traza para cada elemento de la lista.
      
        traz = trace(densit2[i]) #Guardamos el valor de la traza a cada tiempo t.
        
        push!(traza,traz) #Se regresa la lista de tiempo, la lista con la traza al tiempo t, y las matrices GUE 
        #(hamiltoniano libre del ambiente) y la matriz GOE(hamiltoniano de interacción).
    
    end

    return t, traza #Se regresa la lista de tiempo, la lista con la traza al tiempo t.

end

export H3

"""H3(a,n,g,G)

H3 es una función que calcula la pureza de estados generados aleatoriamente con una 
distribución gaussiana a distintos tiempos, considerando dos ambientes (1 y 2) y un potencial de interacción 
sistema-ambiente 1 y otro sistema-ambiente 2, es decir un hamiltoniano de interacción factorizable.
Toma como entradas tres parámetros (a,g,G) y la dimensión (n) del hamiltoniano del ambiente libre de interacciones"""

function H3(a,n,g,G)
    
    Sz = [1 0;0 -1] #El operador matriz de pauli Sigma Z.
    
    goe = herramienta.GOE(2*n) #Definimos el hamiltoniano de interacción entre el sistema y el ambiente_1.
    
    gue = herramienta.GUE(n) #Definimos el hamiltoniano del ambiente libre de interacciones.
    
    t = 0:0.1:100 #Se crea una lista de tiempo.
    
    H = (a/2)*kron(Sz,eye(n),eye(n))+kron(eye(2),gue,eye(n))+g*kron(goe,eye(n))+kron(eye(2),eye(n),herramienta.GOE(n))+G*kron(eye(2),herramienta.GOE(n),herramienta.GOE(n))
    
    V = herramienta.random_state(n*n) #Se crea un estado aleatorio de dimensión nxn.
    
    rho = [] #Se crea una lista para guardar las matrices de densidad en este for se evoluciona el estado 
    #dado en V usando el operador de evolución temporal.
    
    for i in 1:length(t)
        
        psi = (expm(-im*H*t[i]))*(herramienta.projector(V))*(expm(im*H*t[i])) #Se calcula la matriz de densidad al tiempo t del estado V.
        
         psi = psi/norm(psi) #Se normaliza el estado.
   
        push!(rho,psi) #Se guarda en una lista la evolución de la matriz de densidad a distintos tiempos.
    
    end

    partial = [] #Creamos una lista para guardar la traza parcial.
    
    for i in 1:length(rho)
        
        trace = herramienta.partial_trace_pure_bipartite_mat(rho[i],2,1) #En este for se calcula la traza 
        #parcial para la matriz de densidad a cada tiempo.

        push!(partial,trace) #Se guarda la traza parcial en una lista.
    
    end

    densit2 = [] #Se crea una lista para guardar la matriz al cuadrado.

    for i in 1:length(rho) #En este for se calcula el cuadrado de la traza parcial a cada tiempo t de la lista.
        
        densit = (partial[i])*(partial[i]) #Se eleva al cuadrado.
        
        push!(densit2,densit) #Se guarda la matriz obtenida.
    
    end
    
    traza = [] #Finalmente se calcula la traza de las matrices obtenidas en el paso anterior.
    
    for i in 1:length(densit2) #Se calcula la traza para cada elemento de la lista.
        
        traz = trace(densit2[i]) #Guardamos el valor de la traza a cada tiempo t.
        
        push!(traza,traz) #Se regresa la lista de tiempo, la lista con la traza al tiempo t, y las matrices
        #GUE(hamiltoniano libre  del ambiente) y la matriz GOE(hamiltoniano de interacción).
        
    end
    
    return t, traza #Se regresa la lista de tiempo, la lista con la traza al tiempo t.

end

export H4

"""H4(a,n,g,G)

H4 es una función que calcula la pureza de estados generados aleatoriamente con una 
distribución gaussiana a distintos tiempos, considerando dos ambientes (1 y 2) y un potencial de interacción 
sistema-ambiente 1 y otro sistema-ambiente 2, es decir un hamiltoniano de interacción factorizable. En este 
hamiltoniano las dimensiones de los ambientes 1 y 2 son distintas. Toma como entradas tres parámetros (a,g,G)
y la dimensión (n) del hamiltoniano del ambiente libre de interacciones"""


function H4(a,n,m,g,G)
    
    Sz = [1 0;0 -1] #Se crea el operador matriz de pauli Sigma Z.
    
    goe = herramienta.GOE(2*n) #Definimos el hamiltoniano de interacción entre el sistema y el ambiente_1.
    
    gue = herramienta.GUE(n) #Definimos el hamiltoniano del ambiente libre de interacciones.
    
    t = 0:0.1:100 #Se crea una lista de tiempo.
    
    H = (a/2)*kron(Sz,eye(n),eye(m))+kron(eye(2),gue,eye(m))+g*kron(goe,eye(m))+kron(eye(2),eye(n),herramienta.GOE(m))+G*kron(eye(2),herramienta.GOE(n),herramienta.GOE(m))
    
    V = herramienta.random_state(m*n) #Se crea un estado aleatorio de dimensión mxn.
    
    rho = [] #Se crea una lista para guardar las matrices de densidad.
    

    for i in 1:length(t) #En este for se evoluciona el estado dado en V usando el operador de evolución temporal.
        
        psi=(expm(-im*H*t[i]))*(herramienta.projector(V))*(expm(im*H*t[i])) #Se calcula la matriz de densidad 
        #al tiempo t del estado V.
        
       
        
        push!(rho,psi) #Se guarda en una lista la evolución de la matriz de densidad a distintos tiempos.
    
    end

    partial = [] #Creamos una lista para guardar la traza parcial.
    
    for i in 1:length(rho) #En este for se calcula la traza parcial para la matriz de densidad
        #a cada tiempo.

        trace = herramienta.partial_trace_pure_bipartite_mat(rho[i],2,1)

        push!(partial,trace) #Se guarda la traza parcial en una lista.

    end

    densit2 = [] #Se crea una lista para guardar la matriz al cuadrado.

    for i in 1:length(rho) #En este for se calcula el cuadrado de la traza parcial a cada tiempo t de la lista.

        densit = (partial[i])*(partial[i]) #Se eleva al cuadrado.
    
        push!(densit2,densit) #Se guarda la matriz obtenida.

    end

    traza = [] #Finalmente se calcula la traza de las matrices obtenidas en el paso anterior.

    for i in 1:length(densit2) #Se calcula la traza para cada elemento de la lista.

        traz = trace(densit2[i]) #Guardamos el valor de la traza a cada tiempo t.
    
        push!(traza,traz) #Se regresa la lista de tiempo, la lista con la traza al tiempo t, y las 
        #matrices GUE(hamiltoniano libre del ambiente) y la matriz GOE(hamiltoniano de interacción).

    end

    return t, traza #Se regresa la lista de tiempo, la lista con la traza al tiempo t.

end

end




