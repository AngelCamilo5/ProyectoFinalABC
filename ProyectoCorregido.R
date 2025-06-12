tabla = data.frame(
  Variante = character(),
  Paises = character()
)
str(tabla)
obs = list("Omicron", "Estados Unidos, Alemania, India, Turquía, Rusia, Reino Unido,etc"); tabla[1,] = obs
obs = list("Alpha","Reino Unido, Estados Unidos, Alemania, Dinamarca y Suecia"); tabla[2,] = obs
obs = list("Beta","Sudáfrica, Filipinas, Suecia, Alemania y Estados Unidos"); tabla[3,] = obs
obs = list("Gamma","Brasil, Estados Unidos, Chile, Argentina y España"); tabla[4,] = obs
obs = list("Epsilon","Aruba, México, Estados Unidos, Canadá y Argentina"); tabla[5,] = obs
obs = list("Delta","India, Turquía, Estados Unidos, Alemania y Reino Unido"); tabla[6,] = obs
obs = list("Iota","Se encontró por primera vez en Estados Unidos"); tabla[7,] = obs
obs = list("Kappa","Se encontró por primera vez en la India"); tabla[8,] = obs
obs = list("Zeta","Se encontró por primera vez en Brasil"); tabla[9,] = obs
obs = list("Mu","Se encontró por primera vez en Colombia"); tabla[10,] = obs
obs = list("Eta","Se encontró por primera vez en Nigeria y Reino Unido");tabla[11,] = obs
tabla


library(seqinr)
fRef = read.fasta("sequence_ref.txt")
fSouthAfrica1 = read.fasta("Omicron_1_SouthAfrica.fasta")
fFrance1 = read.fasta("Omicron_1_France.fasta")
fIndia1 = read.fasta("Omicron_1_India.fasta")
fMexico1 = read.fasta("Omicron_1_Mexico.fasta")
fUSA1 = read.fasta("Omicron_1_USA.fasta")

library(seqinr)
fSouthAfrica2 = read.fasta("Omicron_2_SouthAfrica.fasta")
fFrance2 = read.fasta("Omicron_2_France.fasta")
fIndia2 = read.fasta("Omicron_2_India.fasta")
fMexico2 = read.fasta("Omicron_2_Mexico.fasta")
fUSA2 = read.fasta("Omicron_2_USA.fasta")

df_Mexico_1 = data.frame (Mutation = character(), 
                          CambioCodon = character(), 
                          CambioAmino = character(),
                          gen = character(),
                          numeroGen = integer()
)
str(df_Mexico_1)

df_Francia_1 = data.frame (Mutation = character(), 
                           CambioCodon = character(), 
                           CambioAmino = character(),
                           gen = character(),
                           numeroGen = integer()
)
str(df_Francia_1)

df_India_1 = data.frame (Mutation = character(), 
                         CambioCodon = character(), 
                         CambioAmino = character(),
                         gen = character(),
                         numeroGen = integer()
)
str(df_India_1)

df_SouthAfrica_1 = data.frame (Mutation = character(), 
                               CambioCodon = character(), 
                               CambioAmino = character(),
                               gen = character(),
                               numeroGen = integer()
)
str(df_SouthAfrica_1)

df_USA_1 = data.frame (Mutation = character(), 
                       CambioCodon = character(), 
                       CambioAmino = character(),
                       gen = character(),
                       numeroGen = integer()
)
str(df_USA_1)


df_Mexico_2 = data.frame (Mutation = character(), 
                          CambioCodon = character(), 
                          CambioAmino = character(),
                          gen = character(),
                          numeroGen = integer()
)
str(df_Mexico_2)

df_Francia_2 = data.frame (Mutation = character(), 
                           CambioCodon = character(), 
                           CambioAmino = character(),
                           gen = character(),
                           numeroGen = integer()
)
str(df_Francia_2)

df_India_2 = data.frame (Mutation = character(), 
                         CambioCodon = character(), 
                         CambioAmino = character(),
                         gen = character(),
                         numeroGen = integer()
)
str(df_India_2)

df_SouthAfrica_2 = data.frame (Mutation = character(), 
                               CambioCodon = character(), 
                               CambioAmino = character(),
                               gen = character(),
                               numeroGen = integer()
)
str(df_SouthAfrica_2)

df_USA_2 = data.frame (Mutation = character(), 
                       CambioCodon = character(), 
                       CambioAmino = character(),
                       gen = character(),
                       numeroGen = integer()
)
str(df_USA_2)


cambio_aminoacido = function(codon) {
  resultado = switch(codon,
                     "GAC" = "D", "GAU" = "D",
                     "GAA" = "E", "GAG" = "E",
                     "CGA" = "R", "CGC" = "R", "CGG" = "R", "CGU" = "R", "AGA" = "R", "AGG" = "R",
                     "AAA" = "K", "AAG" = "K",
                     "AAC" = "N", "AAU" = "N",
                     "CAC" = "H", "CAU" = "H",
                     "CAA" = "Q", "CAG" = "Q",
                     "UCA" = "S", "UCC" = "S", "UCG" = "S", "UCU" = "S", "AGC" = "S", "AGU" = "S",
                     "ACA" = "T", "ACC" = "T", "ACG" = "T", "ACU" = "T",
                     "GCA" = "A", "GCC" = "A", "GCG" = "A", "GCU" = "A",
                     "GGA" = "G", "GGC" = "G", "GGG" = "G", "GGU" = "G",
                     "GUA" = "V", "GUC" = "V", "GUG" = "V", "GUU" = "V",
                     "CCA" = "P", "CCC" = "P", "CCG" = "P", "CCU" = "P",
                     "CUA" = "L", "CUC" = "L", "CUG" = "L", "CUU" = "L", "UUA" = "L", "UUG" = "L",
                     "UUC" = "F", "UUU" = "F",
                     "UAC" = "Y", "UAU" = "Y", 
                     "AUA" = "I", "AUC" = "I", "AUU" = "I",
                     "AUG" = "M",
                     "UGG" = "W",
                     "UGC" = "C", "UGU" = "C",
                     "UAA" = "Terminación", "UAG" = "Terminación", "UGA" = "Terminación")
  return(resultado)
}

CalcularPeso = function(fila, col){
  if (A[col-1]==B[fila-1]) diagonal = m[fila-1,col-1] + 1
  else diagonal = m[fila-1,col-1] - 1
  up = m[fila-1,col] - 2
  left = m[fila,col-1] - 2
  peso = max (diagonal, up, left)
  return (peso)
}


cant = 0
cuenta = 0
genRef = as.vector( fRef[[6]] )
for (i in seq(6, length(fMexico1),12)){
  genMexico1 = as.vector( fMexico1[[i]] )
  genMexico1[which(genMexico1=="n")] = "x"
  dif = which(genRef != genMexico1)
  genRef[which(genRef=="t")] = "u"
  genMexico1[which(genMexico1=="t")] = "u"
  cuenta = cuenta + 1
  if (length(genRef)==length(genMexico1)){
    dif = which(genRef != genMexico1)
    for (j in dif) {
      mutacion = paste(toupper(genRef[j]), "to", toupper(genMexico1[j]), sep="")
      num_codon= ((j-1)%/%3)+1
      pos_en_codon=((j-1)%%3)+1
      if (pos_en_codon==1){
        anterior = paste(genRef[j:(j+2)], collapse = "")
        nuevo = paste(genMexico1[j:(j+2)], collapse = "")
        cambio_en_codon = paste(toupper(anterior), "to", toupper(nuevo), sep="")
      }
      if (pos_en_codon==2){
        anterior = paste(genRef[(j-1):(j+1)], collapse = "")
        nuevo = paste(genMexico1[(j-1):(j+1)], collapse = "")
        cambio_en_codon = paste(toupper(anterior), "to", toupper(nuevo), sep="")
      }
      if (pos_en_codon==3){
        anterior = paste(genRef[(j-2):(j)], collapse = "")
        nuevo = paste(genMexico1[(j-2):(j)], collapse = "")
        cambio_en_codon = paste(toupper(anterior), "to", toupper(nuevo), sep="")
      }
      if (!grepl("x", nuevo)) {
        amin1 = cambio_aminoacido(toupper(anterior))
        amin2 = cambio_aminoacido(toupper(nuevo))
        if (amin1 != amin2) {
          cambio_en_amin = paste(amin1, num_codon, amin2, sep="")
          cant = cant + 1
          obs = list(mutacion,cambio_en_codon,cambio_en_amin,"M", cuenta); df_Mexico_1[cant,] = obs
        }
      }
    }
  }
  else {
    cant = cant + 1
    obs = list("Deleción o inersción","Deleción o inersción","Deleción o inersción","M", cuenta);
    df_Mexico_1[cant,] = obs
  }
}


cant = 0
cuenta = 0
genRef = as.vector( fRef[[6]] )
for (i in seq(6, length(fFrance1),12)){
  genFrance1 = as.vector( fFrance1[[i]] )
  genFrance1[which(genFrance1=="n")] = "x"
  dif = which(genRef != genFrance1)
  genRef[which(genRef=="t")] = "u"
  genFrance1[which(genFrance1=="t")] = "u"
  genFrance1[which(genFrance1=="r")] = "x"
  genFrance1[which(genFrance1=="s")] = "x"
  cuenta = cuenta + 1
  if (length(genRef)==length(genFrance1)){
    dif = which(genRef != genFrance1)
    for (j in dif) {
      mutacion = paste(toupper(genRef[j]), "to", toupper(genFrance1[j]), sep="")
      num_codon= ((j-1)%/%3)+1
      pos_en_codon=((j-1)%%3)+1
      if (pos_en_codon==1){
        anterior = paste(genRef[j:(j+2)], collapse = "")
        nuevo = paste(genFrance1[j:(j+2)], collapse = "")
        cambio_en_codon = paste(toupper(anterior), "to", toupper(nuevo), sep="")
      }
      if (pos_en_codon==2){
        anterior = paste(genRef[(j-1):(j+1)], collapse = "")
        nuevo = paste(genFrance1[(j-1):(j+1)], collapse = "")
        cambio_en_codon = paste(toupper(anterior), "to", toupper(nuevo), sep="")
      }
      if (pos_en_codon==3){
        anterior = paste(genRef[(j-2):(j)], collapse = "")
        nuevo = paste(genFrance1[(j-2):(j)], collapse = "")
        cambio_en_codon = paste(toupper(anterior), "to", toupper(nuevo), sep="")
      }
      if (!grepl("x", nuevo)) {
        amin1 = cambio_aminoacido(toupper(anterior))
        amin2 = cambio_aminoacido(toupper(nuevo))
        if (amin1 != amin2) {
          cambio_en_amin = paste(amin1, num_codon, amin2, sep="")
          cant = cant + 1
          obs = list(mutacion,cambio_en_codon,cambio_en_amin,"M", cuenta); df_Francia_1[cant,] = obs
        }
      }
    }
  }
  else {
    cant = cant + 1
    obs = list("Deleción o inersción","Deleción o inersción","Deleción o inersción","M", cuenta);
    df_Francia_1[cant,] = obs
  }
}


cant = 0
cuenta = 0
genRef = as.vector( fRef[[6]] )
for (i in seq(6, length(fUSA1),12)){
  genUSA1 = as.vector( fUSA1[[i]] )
  genUSA1[which(genUSA1=="n")] = "x"
  dif = which(genRef != genUSA1)
  genRef[which(genRef=="t")] = "u"
  genUSA1[which(genUSA1=="t")] = "u"
  cuenta = cuenta + 1
  if (length(genRef)==length(genUSA1)){
    dif = which(genRef != genUSA1)
    for (j in dif) {
      mutacion = paste(toupper(genRef[j]), "to", toupper(genUSA1[j]), sep="")
      num_codon= ((j-1)%/%3)+1
      pos_en_codon=((j-1)%%3)+1
      if (pos_en_codon==1){
        anterior = paste(genRef[j:(j+2)], collapse = "")
        nuevo = paste(genUSA1[j:(j+2)], collapse = "")
        cambio_en_codon = paste(toupper(anterior), "to", toupper(nuevo), sep="")
      }
      if (pos_en_codon==2){
        anterior = paste(genRef[(j-1):(j+1)], collapse = "")
        nuevo = paste(genUSA1[(j-1):(j+1)], collapse = "")
        cambio_en_codon = paste(toupper(anterior), "to", toupper(nuevo), sep="")
      }
      if (pos_en_codon==3){
        anterior = paste(genRef[(j-2):(j)], collapse = "")
        nuevo = paste(genUSA1[(j-2):(j)], collapse = "")
        cambio_en_codon = paste(toupper(anterior), "to", toupper(nuevo), sep="")
      }
      if (!grepl("x", nuevo)) {
        amin1 = cambio_aminoacido(toupper(anterior))
        amin2 = cambio_aminoacido(toupper(nuevo))
        if (amin1 != amin2) {
          cambio_en_amin = paste(amin1, num_codon, amin2, sep="")
          cant = cant + 1
          obs = list(mutacion,cambio_en_codon,cambio_en_amin,"M", cuenta); df_USA_1[cant,] = obs
        }
      }
    }
  }
  else {
    cant = cant + 1
    obs = list("Deleción o inersción","Deleción o inersción","Deleción o inersción","M", cuenta);
    df_USA_1[cant,] = obs
  }
}


cant = 0
cuenta = 0
genRef = as.vector( fRef[[6]] )
for (i in seq(6, length(fIndia1),12)){
  genIndia1 = as.vector( fIndia1[[i]] )
  genIndia1[which(genIndia1=="n")] = "x"
  dif = which(genRef != genIndia1)
  genRef[which(genRef=="t")] = "u"
  genIndia1[which(genIndia1=="t")] = "u"
  cuenta = cuenta + 1
  if (length(genRef)==length(genIndia1)){
    dif = which(genRef != genIndia1)
    for (j in dif) {
      mutacion = paste(toupper(genRef[j]), "to", toupper(genIndia1[j]), sep="")
      num_codon= ((j-1)%/%3)+1
      pos_en_codon=((j-1)%%3)+1
      if (pos_en_codon==1){
        anterior = paste(genRef[j:(j+2)], collapse = "")
        nuevo = paste(genIndia1[j:(j+2)], collapse = "")
        cambio_en_codon = paste(toupper(anterior), "to", toupper(nuevo), sep="")
      }
      if (pos_en_codon==2){
        anterior = paste(genRef[(j-1):(j+1)], collapse = "")
        nuevo = paste(genIndia1[(j-1):(j+1)], collapse = "")
        cambio_en_codon = paste(toupper(anterior), "to", toupper(nuevo), sep="")
      }
      if (pos_en_codon==3){
        anterior = paste(genRef[(j-2):(j)], collapse = "")
        nuevo = paste(genIndia1[(j-2):(j)], collapse = "")
        cambio_en_codon = paste(toupper(anterior), "to", toupper(nuevo), sep="")
      }
      if (!grepl("x", nuevo)) {
        amin1 = cambio_aminoacido(toupper(anterior))
        amin2 = cambio_aminoacido(toupper(nuevo))
        if (amin1 != amin2) {
          cambio_en_amin = paste(amin1, num_codon, amin2, sep="")
          cant = cant + 1
          obs = list(mutacion,cambio_en_codon,cambio_en_amin,"M", cuenta); df_India_1[cant,] = obs
        }
      }
    }
  }
  else {
    cant = cant + 1
    obs = list("Deleción o inersción","Deleción o inersción","Deleción o inersción","M", cuenta);
    df_India_1[cant,] = obs
  }
}


cant = 0
cuenta = 0
genRef = as.vector( fRef[[6]] )
for (i in seq(6, length(fSouthAfrica1),12)){
  genSouthAfrica1 = as.vector( fSouthAfrica1[[i]] )
  genSouthAfrica1[which(genSouthAfrica1=="n")] = "x"
  dif = which(genRef != genSouthAfrica1)
  genRef[which(genRef=="t")] = "u"
  genSouthAfrica1[which(genSouthAfrica1=="t")] = "u"
  cuenta = cuenta + 1
  if (length(genRef)==length(genSouthAfrica1)){
    dif = which(genRef != genSouthAfrica1)
    for (j in dif) {
      mutacion = paste(toupper(genRef[j]), "to", toupper(genSouthAfrica1[j]), sep="")
      num_codon= ((j-1)%/%3)+1
      pos_en_codon=((j-1)%%3)+1
      if (pos_en_codon==1){
        anterior = paste(genRef[j:(j+2)], collapse = "")
        nuevo = paste(genSouthAfrica1[j:(j+2)], collapse = "")
        cambio_en_codon = paste(toupper(anterior), "to", toupper(nuevo), sep="")
      }
      if (pos_en_codon==2){
        anterior = paste(genRef[(j-1):(j+1)], collapse = "")
        nuevo = paste(genSouthAfrica1[(j-1):(j+1)], collapse = "")
        cambio_en_codon = paste(toupper(anterior), "to", toupper(nuevo), sep="")
      }
      if (pos_en_codon==3){
        anterior = paste(genRef[(j-2):(j)], collapse = "")
        nuevo = paste(genSouthAfrica1[(j-2):(j)], collapse = "")
        cambio_en_codon = paste(toupper(anterior), "to", toupper(nuevo), sep="")
      }
      if (!grepl("x", nuevo)) {
        amin1 = cambio_aminoacido(toupper(anterior))
        amin2 = cambio_aminoacido(toupper(nuevo))
        if (amin1 != amin2) {
          cambio_en_amin = paste(amin1, num_codon, amin2, sep="")
          cant = cant + 1
          obs = list(mutacion,cambio_en_codon,cambio_en_amin,"M", cuenta); df_SouthAfrica_1[cant,] = obs
        }
      }
    }
  }
  else {
    cant = cant + 1
    obs = list("Deleción o inersción","Deleción o inersción","Deleción o inersción","M", cuenta);
    df_SouthAfrica_1[cant,] = obs
  }
}


cant = 0
cuenta = 0
genRef = as.vector( fRef[[6]] )
for (i in seq(6, length(fMexico2),12)){
  genMexico2 = as.vector( fMexico2[[i]] )
  genMexico2[which(genMexico2=="n")] = "x"
  dif = which(genRef != genMexico2)
  genRef[which(genRef=="t")] = "u"
  genMexico2[which(genMexico2=="t")] = "u"
  cuenta = cuenta + 1
  if (length(genRef)==length(genMexico2)){
    dif = which(genRef != genMexico2)
    for (j in dif) {
      mutacion = paste(toupper(genRef[j]), "to", toupper(genMexico2[j]), sep="")
      num_codon= ((j-1)%/%3)+1
      pos_en_codon=((j-1)%%3)+1
      if (pos_en_codon==1){
        anterior = paste(genRef[j:(j+2)], collapse = "")
        nuevo = paste(genMexico2[j:(j+2)], collapse = "")
        cambio_en_codon = paste(toupper(anterior), "to", toupper(nuevo), sep="")
      }
      if (pos_en_codon==2){
        anterior = paste(genRef[(j-1):(j+1)], collapse = "")
        nuevo = paste(genMexico2[(j-1):(j+1)], collapse = "")
        cambio_en_codon = paste(toupper(anterior), "to", toupper(nuevo), sep="")
      }
      if (pos_en_codon==3){
        anterior = paste(genRef[(j-2):(j)], collapse = "")
        nuevo = paste(genMexico2[(j-2):(j)], collapse = "")
        cambio_en_codon = paste(toupper(anterior), "to", toupper(nuevo), sep="")
      }
      if (!grepl("x", nuevo)) {
        amin1 = cambio_aminoacido(toupper(anterior))
        amin2 = cambio_aminoacido(toupper(nuevo))
        if (amin1 != amin2) {
          cambio_en_amin = paste(amin1, num_codon, amin2, sep="")
          cant = cant + 1
          obs = list(mutacion,cambio_en_codon,cambio_en_amin,"M", cuenta); df_Mexico_2[cant,] = obs
        }
      }
    }
  }
  else {
    cant = cant + 1
    obs = list("Deleción o inersción","Deleción o inersción","Deleción o inersción","M", cuenta);
    df_Mexico_2[cant,] = obs
  }
}


cant = 0
cuenta = 0
genRef = as.vector( fRef[[6]] )
for (i in seq(6, length(fFrance2),12)){
  genFrance2 = as.vector( fFrance2[[i]] )
  genFrance2[which(genFrance2=="n")] = "x"
  dif = which(genRef != genFrance2)
  genRef[which(genRef=="t")] = "u"
  genFrance2[which(genFrance2=="t")] = "u"
  genFrance2[which(genFrance2=="r")] = "x"
  genFrance2[which(genFrance2=="s")] = "x"
  cuenta = cuenta + 1
  if (length(genRef)==length(genFrance2)){
    dif = which(genRef != genFrance2)
    for (j in dif) {
      mutacion = paste(toupper(genRef[j]), "to", toupper(genFrance2[j]), sep="")
      num_codon= ((j-1)%/%3)+1
      pos_en_codon=((j-1)%%3)+1
      if (pos_en_codon==1){
        anterior = paste(genRef[j:(j+2)], collapse = "")
        nuevo = paste(genFrance2[j:(j+2)], collapse = "")
        cambio_en_codon = paste(toupper(anterior), "to", toupper(nuevo), sep="")
      }
      if (pos_en_codon==2){
        anterior = paste(genRef[(j-1):(j+1)], collapse = "")
        nuevo = paste(genFrance2[(j-1):(j+1)], collapse = "")
        cambio_en_codon = paste(toupper(anterior), "to", toupper(nuevo), sep="")
      }
      if (pos_en_codon==3){
        anterior = paste(genRef[(j-2):(j)], collapse = "")
        nuevo = paste(genFrance2[(j-2):(j)], collapse = "")
        cambio_en_codon = paste(toupper(anterior), "to", toupper(nuevo), sep="")
      }
      if (!grepl("x", nuevo)) {
        amin1 = cambio_aminoacido(toupper(anterior))
        amin2 = cambio_aminoacido(toupper(nuevo))
        if (amin1 != amin2) {
          cambio_en_amin = paste(amin1, num_codon, amin2, sep="")
          cant = cant + 1
          obs = list(mutacion,cambio_en_codon,cambio_en_amin,"M", cuenta); df_Francia_2[cant,] = obs
        }
      }
    }
  }
  else {
    cant = cant + 1
    obs = list("Deleción o inersción","Deleción o inersción","Deleción o inersción","M", cuenta);
    df_Francia_2[cant,] = obs
  }
}


cant = 0
cuenta = 0
genRef = as.vector( fRef[[6]] )
for (i in seq(6, length(fUSA2),12)){
  genUSA2 = as.vector( fUSA2[[i]] )
  genUSA2[which(genUSA2=="n")] = "x"
  genUSA2[which(genUSA2=="r")] = "x"
  genUSA2[which(genUSA2=="y")] = "x"
  dif = which(genRef != genUSA2)
  genRef[which(genRef=="t")] = "u"
  genUSA2[which(genUSA2=="t")] = "u"
  cuenta = cuenta + 1
  if (length(genRef)==length(genUSA2)){
    dif = which(genRef != genUSA2)
    for (j in dif) {
      mutacion = paste(toupper(genRef[j]), "to", toupper(genUSA2[j]), sep="")
      num_codon= ((j-1)%/%3)+1
      pos_en_codon=((j-1)%%3)+1
      if (pos_en_codon==1){
        anterior = paste(genRef[j:(j+2)], collapse = "")
        nuevo = paste(genUSA2[j:(j+2)], collapse = "")
        cambio_en_codon = paste(toupper(anterior), "to", toupper(nuevo), sep="")
      }
      if (pos_en_codon==2){
        anterior = paste(genRef[(j-1):(j+1)], collapse = "")
        nuevo = paste(genUSA2[(j-1):(j+1)], collapse = "")
        cambio_en_codon = paste(toupper(anterior), "to", toupper(nuevo), sep="")
      }
      if (pos_en_codon==3){
        anterior = paste(genRef[(j-2):(j)], collapse = "")
        nuevo = paste(genUSA2[(j-2):(j)], collapse = "")
        cambio_en_codon = paste(toupper(anterior), "to", toupper(nuevo), sep="")
      }
      if (!grepl("x", nuevo)) {
        amin1 = cambio_aminoacido(toupper(anterior))
        amin2 = cambio_aminoacido(toupper(nuevo))
        if (amin1 != amin2) {
          cambio_en_amin = paste(amin1, num_codon, amin2, sep="")
          cant = cant + 1
          obs = list(mutacion,cambio_en_codon,cambio_en_amin,"M", cuenta); df_USA_2[cant,] = obs
        }
      }
    }
  }
  else {
    cant = cant + 1
    obs = list("Deleción o inersción","Deleción o inersción","Deleción o inersción","M", cuenta);
    df_USA_2[cant,] = obs
  }
}


cant = 0
cuenta = 0
genRef = as.vector( fRef[[6]] )
for (i in seq(6, length(fIndia2),12)){
  genIndia2 = as.vector( fIndia2[[i]] )
  genIndia2[which(genIndia2=="n")] = "x"
  genIndia2[which(genIndia2=="y")] = "x"
  genIndia2[which(genIndia2=="r")] = "x"
  genIndia2[which(genIndia2=="s")] = "x"
  dif = which(genRef != genIndia2)
  genRef[which(genRef=="t")] = "u"
  genIndia2[which(genIndia2=="t")] = "u"
  cuenta = cuenta + 1
  if (length(genRef)==length(genIndia2)){
    dif = which(genRef != genIndia2)
    for (j in dif) {
      mutacion = paste(toupper(genRef[j]), "to", toupper(genIndia2[j]), sep="")
      num_codon= ((j-1)%/%3)+1
      pos_en_codon=((j-1)%%3)+1
      if (pos_en_codon==1){
        anterior = paste(genRef[j:(j+2)], collapse = "")
        nuevo = paste(genIndia2[j:(j+2)], collapse = "")
        cambio_en_codon = paste(toupper(anterior), "to", toupper(nuevo), sep="")
      }
      if (pos_en_codon==2){
        anterior = paste(genRef[(j-1):(j+1)], collapse = "")
        nuevo = paste(genIndia2[(j-1):(j+1)], collapse = "")
        cambio_en_codon = paste(toupper(anterior), "to", toupper(nuevo), sep="")
      }
      if (pos_en_codon==3){
        anterior = paste(genRef[(j-2):(j)], collapse = "")
        nuevo = paste(genIndia2[(j-2):(j)], collapse = "")
        cambio_en_codon = paste(toupper(anterior), "to", toupper(nuevo), sep="")
      }
      if (!grepl("x", nuevo)) {
        amin1 = cambio_aminoacido(toupper(anterior))
        amin2 = cambio_aminoacido(toupper(nuevo))
        if (amin1 != amin2) {
          cambio_en_amin = paste(amin1, num_codon, amin2, sep="")
          cant = cant + 1
          obs = list(mutacion,cambio_en_codon,cambio_en_amin,"M", cuenta); df_India_2[cant,] = obs
        }
      }
    }
  }
  else {
    cant = cant + 1
    obs = list("Deleción o inersción","Deleción o inersción","Deleción o inersción","M", cuenta);
    df_India_2[cant,] = obs
  }
}


cant = 0
cuenta = 0
genRef = as.vector( fRef[[6]] )
for (i in seq(6, length(fSouthAfrica2),12)){
  genSouthAfrica2 = as.vector( fSouthAfrica2[[i]] )
  genSouthAfrica2[which(genSouthAfrica2=="n")] = "x"
  dif = which(genRef != genSouthAfrica2)
  genRef[which(genRef=="t")] = "u"
  genSouthAfrica2[which(genSouthAfrica2=="t")] = "u"
  cuenta = cuenta + 1
  if (length(genRef)==length(genSouthAfrica2)){
    dif = which(genRef != genSouthAfrica2)
    for (j in dif) {
      mutacion = paste(toupper(genRef[j]), "to", toupper(genSouthAfrica2[j]), sep="")
      num_codon= ((j-1)%/%3)+1
      pos_en_codon=((j-1)%%3)+1
      if (pos_en_codon==1){
        anterior = paste(genRef[j:(j+2)], collapse = "")
        nuevo = paste(genSouthAfrica2[j:(j+2)], collapse = "")
        cambio_en_codon = paste(toupper(anterior), "to", toupper(nuevo), sep="")
      }
      if (pos_en_codon==2){
        anterior = paste(genRef[(j-1):(j+1)], collapse = "")
        nuevo = paste(genSouthAfrica2[(j-1):(j+1)], collapse = "")
        cambio_en_codon = paste(toupper(anterior), "to", toupper(nuevo), sep="")
      }
      if (pos_en_codon==3){
        anterior = paste(genRef[(j-2):(j)], collapse = "")
        nuevo = paste(genSouthAfrica2[(j-2):(j)], collapse = "")
        cambio_en_codon = paste(toupper(anterior), "to", toupper(nuevo), sep="")
      }
      if (!grepl("x", nuevo)) {
        amin1 = cambio_aminoacido(toupper(anterior))
        amin2 = cambio_aminoacido(toupper(nuevo))
        if (amin1 != amin2) {
          cambio_en_amin = paste(amin1, num_codon, amin2, sep="")
          cant = cant + 1
          obs = list(mutacion,cambio_en_codon,cambio_en_amin,"M", cuenta); df_SouthAfrica_2[cant,] = obs
        }
      }
    }
  }
  else {
    cant = cant + 1
    obs = list("Deleción o inersción","Deleción o inersción","Deleción o inersción","M", cuenta);
    df_SouthAfrica_2[cant,] = obs
  }
}



library(dplyr)
library(ggplot2)
df_Mexico_Amino1 = filter(
  dplyr::summarize(
    select(
      group_by(df_Mexico_1,CambioAmino),
      Mutation:gen
    ),
    Mutation = first(Mutation),
    CambioCodon = first(CambioCodon),
    gen = first(gen),
    Cuenta = n()
  ),
  Cuenta>=1
)

df_Mexico_Amino1 = df_Mexico_Amino1[order(-df_Mexico_Amino1$Cuenta), ]
df_Mexico_Amino1


pM1 = ggplot(df_Mexico_Amino1)
pM1 = pM1 + aes(x=CambioAmino, y=Cuenta, fill=CambioAmino, label=Cuenta)
pM1 = pM1 + ggtitle("Cambio de Aminoácidos BA.1")
pM1 = pM1 + labs(x="Amino", y="Frecuencia", fill="Frecuencia")
pM1 = pM1 + geom_bar(stat = "identity")
pM1 = pM1 + geom_text(stat = "identity", vjust=1.5)
pM1


library(dplyr)
df_Francia_Amino1 = filter(
  dplyr::summarize(
    select(
      group_by(df_Francia_1,CambioAmino),
      Mutation:gen
    ),
    Mutation = first(Mutation),
    CambioCodon = first(CambioCodon),
    gen = first(gen),
    Cuenta = n()
  ),
  Cuenta>=1
)

df_Francia_Amino1 = df_Francia_Amino1[order(-df_Francia_Amino1$Cuenta), ]
df_Francia_Amino1


pF1 = ggplot(df_Francia_Amino1)
pF1 = pF1 + aes(x=CambioAmino, y=Cuenta, fill=CambioAmino, label=Cuenta)
pF1 = pF1 + ggtitle("Cambio de Aminoácidos BA.1")
pF1 = pF1 + labs(x="Amino", y="Frecuencia", fill="Frecuencia")
pF1 = pF1 + geom_bar(stat = "identity")
pF1 = pF1 + geom_text(stat = "identity", vjust=1.5)
pF1


library(dplyr)
df_USA_Amino1 = filter(
  dplyr::summarize(
    select(
      group_by(df_USA_1,CambioAmino),
      Mutation:gen
    ),
    Mutation = first(Mutation),
    CambioCodon = first(CambioCodon),
    gen = first(gen),
    Cuenta = n()
  ),
  Cuenta>=1
)

df_USA_Amino1 = df_USA_Amino1[order(-df_USA_Amino1$Cuenta), ]
df_USA_Amino1

pUSA1 = ggplot(df_USA_Amino1)
pUSA1 = pUSA1 + aes(x=CambioAmino, y=Cuenta, fill=CambioAmino, label=Cuenta)
pUSA1 = pUSA1 + ggtitle("Cambio de Aminoácidos BA.1")
pUSA1 = pUSA1 + labs(x="Amino", y="Frecuencia", fill="Frecuencia")
pUSA1 = pUSA1 + geom_bar(stat = "identity")
pUSA1 = pUSA1 + geom_text(stat = "identity", vjust=1.5)
pUSA1


library(dplyr)
df_India_Amino1 = filter(
  dplyr::summarize(
    select(
      group_by(df_India_1,CambioAmino),
      Mutation:gen
    ),
    Mutation = first(Mutation),
    CambioCodon = first(CambioCodon),
    gen = first(gen),
    Cuenta = n()
  ),
  Cuenta>=1
)

df_India_Amino1 = df_India_Amino1[order(-df_India_Amino1$Cuenta), ]
df_India_Amino1


pI1 = ggplot(df_India_Amino1)
pI1 = pI1 + aes(x=CambioAmino, y=Cuenta, fill=CambioAmino, label=Cuenta)
pI1 = pI1 + ggtitle("Cambio de Aminoácidos BA.1")
pI1 = pI1 + labs(x="Amino", y="Frecuencia", fill="Frecuencia")
pI1 = pI1 + geom_bar(stat = "identity")
pI1 = pI1 + geom_text(stat = "identity", vjust=1.5)
pI1

library(dplyr)
df_SouthAfrica_Amino1 = filter(
  dplyr::summarize(
    select(
      group_by(df_SouthAfrica_1,CambioAmino),
      Mutation:gen
    ),
    Mutation = first(Mutation),
    CambioCodon = first(CambioCodon),
    gen = first(gen),
    Cuenta = n()
  ),
  Cuenta>=1
)

df_SouthAfrica_Amino1 = df_SouthAfrica_Amino1[order(-df_SouthAfrica_Amino1$Cuenta), ]
df_SouthAfrica_Amino1

pSA1 = ggplot(df_SouthAfrica_Amino1)
pSA1 = pSA1 + aes(x=CambioAmino, y=Cuenta, fill=CambioAmino, label=Cuenta)
pSA1 = pSA1 + ggtitle("Cambio de Aminoácidos BA.1")
pSA1 = pSA1 + labs(x="Amino", y="Frecuencia", fill="Frecuencia")
pSA1 = pSA1 + geom_bar(stat = "identity")
pSA1 = pSA1 + geom_text(stat = "identity", vjust=1.5)
pSA1


tabla_combinada_BA1 <- rbind(df_Mexico_1, df_Francia_1, df_India_1, df_SouthAfrica_1, df_USA_1)

library(dplyr)
tabla_final_BA1 = filter(
  summarise(
    select(
      group_by(tabla_combinada_BA1, CambioAmino),
      Mutation:gen
    ),
    Mutation = first(Mutation),
    CambioAmino = first(CambioAmino),
    gen = first(gen),
    Cuenta = n()
  ),
  Cuenta >= 1 
)

tabla_final_BA1 = tabla_final_BA1[order(-tabla_final_BA1$Cuenta), ]
tabla_final_BA1

pBA1 = ggplot(tabla_final_BA1)
pBA1 = pBA1 + aes(x=CambioAmino, y=Cuenta, fill=CambioAmino, label=Cuenta)
pBA1 = pBA1 + ggtitle("Cambio de Aminoácidos BA.1")
pBA1 = pBA1 + labs(x="Amino", y="Frecuencia", fill="Frecuencia")
pBA1 = pBA1 + geom_bar(stat = "identity")
pBA1 = pBA1 + geom_text(stat = "identity", vjust=1.5)
pBA1


library(dplyr)
df_Mexico_Amino2 = filter(
  dplyr::summarize(
    select(
      group_by(df_Mexico_2,CambioAmino),
      Mutation:gen
    ),
    Mutation = first(Mutation),
    CambioCodon = first(CambioCodon),
    gen = first(gen),
    Cuenta = n()
  ),
  Cuenta>=1
)

df_Mexico_Amino2 = df_Mexico_Amino2[order(-df_Mexico_Amino2$Cuenta), ]
df_Mexico_Amino2

pM2 = ggplot(df_Mexico_Amino2)
pM2 = pM2 + aes(x=CambioAmino, y=Cuenta, fill=CambioAmino, label=Cuenta)
pM2 = pM2 + ggtitle("Cambio de Aminoácidos BA.2")
pM2 = pM2 + labs(x="Amino", y="Frecuencia", fill="Frecuencia")
pM2 = pM2 + geom_bar(stat = "identity")
pM2 = pM2 + geom_text(stat = "identity", vjust=1.5)
pM2


library(dplyr)
df_Francia_Amino2 = filter(
  dplyr::summarize(
    select(
      group_by(df_Francia_2,CambioAmino),
      Mutation:gen
    ),
    Mutation = first(Mutation),
    CambioCodon = first(CambioCodon),
    gen = first(gen),
    Cuenta = n()
  ),
  Cuenta>=1
)

df_Francia_Amino2 = df_Francia_Amino2[order(-df_Francia_Amino2$Cuenta), ]
df_Francia_Amino2


pF2 = ggplot(df_Francia_Amino2)
pF2 = pF2 + aes(x=CambioAmino, y=Cuenta, fill=CambioAmino, label=Cuenta)
pF2 = pF2 + ggtitle("Cambio de Aminoácidos BA.2")
pF2 = pF2 + labs(x="Amino", y="Frecuencia", fill="Frecuencia")
pF2 = pF2 + geom_bar(stat = "identity")
pF2 = pF2 + geom_text(stat = "identity", vjust=1.5)
pF2

library(dplyr)
df_USA_Amino2 = filter(
  dplyr::summarize(
    select(
      group_by(df_USA_2,CambioAmino),
      Mutation:gen
    ),
    Mutation = first(Mutation),
    CambioCodon = first(CambioCodon),
    gen = first(gen),
    Cuenta = n()
  ),
  Cuenta>=1
)

df_USA_Amino2 = df_USA_Amino2[order(-df_USA_Amino2$Cuenta), ]
df_USA_Amino2

pUSA2 = ggplot(df_USA_Amino2)
pUSA2 = pUSA2 + aes(x=CambioAmino, y=Cuenta, fill=CambioAmino, label=Cuenta)
pUSA2 = pUSA2 + ggtitle("Cambio de Aminoácidos BA.2")
pUSA2 = pUSA2 + labs(x="Amino", y="Frecuencia", fill="Frecuencia")
pUSA2 = pUSA2 + geom_bar(stat = "identity")
pUSA2 = pUSA2 + geom_text(stat = "identity", vjust=1.5)
pUSA2


library(dplyr)
df_India_Amino2 = filter(
  dplyr::summarize(
    select(
      group_by(df_India_2,CambioAmino),
      Mutation:gen
    ),
    Mutation = first(Mutation),
    CambioCodon = first(CambioCodon),
    gen = first(gen),
    Cuenta = n()
  ),
  Cuenta>=1
)

df_India_Amino2 = df_India_Amino2[order(-df_India_Amino2$Cuenta), ]
df_India_Amino2

pI2 = ggplot(df_India_Amino2)
pI2 = pI2 + aes(x=CambioAmino, y=Cuenta, fill=CambioAmino, label=Cuenta)
pI2 = pI2 + ggtitle("Cambio de Aminoácidos BA.2")
pI2 = pI2 + labs(x="Amino", y="Frecuencia", fill="Frecuencia")
pI2 = pI2 + geom_bar(stat = "identity")
pI2 = pI2 + geom_text(stat = "identity", vjust=1.5)
pI2


library(dplyr)
df_SouthAfrica_Amino2 = filter(
  dplyr::summarize(
    select(
      group_by(df_SouthAfrica_2,CambioAmino),
      Mutation:gen
    ),
    Mutation = first(Mutation),
    CambioCodon = first(CambioCodon),
    gen = first(gen),
    Cuenta = n()
  ),
  Cuenta>=1
)

df_SouthAfrica_Amino2 = df_SouthAfrica_Amino2[order(-df_SouthAfrica_Amino2$Cuenta), ]
df_SouthAfrica_Amino2

pSA2 = ggplot(df_SouthAfrica_Amino2)
pSA2 = pSA2 + aes(x=CambioAmino, y=Cuenta, fill=CambioAmino, label=Cuenta)
pSA2 = pSA2 + ggtitle("Cambio de Aminoácidos BA.2")
pSA2 = pSA2 + labs(x="Amino", y="Frecuencia", fill="Frecuencia")
pSA2 = pSA2 + geom_bar(stat = "identity")
pSA2 = pSA2 + geom_text(stat = "identity", vjust=1.5)
pSA2


tabla_combinada_BA2 <- rbind(df_Mexico_2, df_Francia_2, df_India_2, df_SouthAfrica_2, df_USA_2)

library(dplyr)
tabla_final_BA2 = filter(
  summarise(
    select(
      group_by(tabla_combinada_BA2, CambioAmino),
      Mutation:gen
    ),
    Mutation = first(Mutation),
    CambioAmino = first(CambioAmino),
    gen = first(gen),
    Cuenta = n()
  ),
  Cuenta >= 1 
)

tabla_final_BA2 = tabla_final_BA2[order(-tabla_final_BA2$Cuenta), ]
tabla_final_BA2


pBA2 = ggplot(tabla_final_BA2)
pBA2 = pBA2 + aes(x=CambioAmino, y=Cuenta, fill=CambioAmino, label=Cuenta)
pBA2 = pBA2 + ggtitle("Cambio de Aminoácidos BA.2")
pBA2 = pBA2 + labs(x="Amino", y="Frecuencia", fill="Frecuencia")
pBA2 = pBA2 + geom_bar(stat = "identity")
pBA2 = pBA2 + geom_text(stat = "identity", vjust=1.5)
pBA2


tabla_combinada_linajes <- rbind(tabla_final_BA1, tabla_final_BA2)

library(dplyr)
tabla_final_linajes = filter(
  summarise(
    select(
      group_by(tabla_combinada_linajes, CambioAmino),
      Mutation:gen
    ),
    Mutation = first(Mutation),
    CambioAmino = first(CambioAmino),
    gen = first(gen),
    Cuenta = n()
  ),
  Cuenta >= 1 & Cuenta < 3
)

tabla_final_linajes = tabla_final_linajes[order(-tabla_final_linajes$Cuenta), ]
tabla_final_linajes


pFinal = ggplot(tabla_final_linajes)
pFinal = pFinal + aes(x=CambioAmino, y=Cuenta, fill=CambioAmino, label=Cuenta)
pFinal = pFinal + ggtitle("Cambio de Aminoácidos Ambos")
pFinal = pFinal + labs(x="Amino", y="Frecuencia", fill="Frecuencia")
pFinal = pFinal + geom_bar(stat = "identity")
pFinal = pFinal + geom_text(stat = "identity", vjust=1.5)
pFinal
