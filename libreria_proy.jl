
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

end




