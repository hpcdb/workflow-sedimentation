#Caminho até o jar do RDI
export RDI_DIR='/home/thaylon/Commands/RDI-1.0-jar-with-dependencies.jar'

#Criando o comando "RDI"
alias RDI="java -jar \"$RDI_DIR\""

#Caminho até diretório bin/ do Fastbit
export FASTBIT_BIN='/home/thaylon/Softwares/fastbit-2.0.3/bin/'

#As colunas presente no arquivo CSV
#Existem 3 tipos de colunas, sendo eles
#Arquivo -> file
#Texto -> text
#Número inteiro -> number
#Número decimal -> numeric:<número de casas depois da vírgula>, para o fastbit essa precisão não importa
export COLUMNS='[u:numeric:5,v:numeric:5,w:numeric:5,p:numeric:5,s:numeric:5,d:numeric:5,vtkValidPointMask:numeric,arc_length:numeric:5,Points0:numeric:5,Points1:numeric:5,Points2:numeric:5]'

#Indexação por regex, no caso todos arquivos .csv
export RDI_CMD="$PWD \"*.csv\" $COLUMNS -bin=\"$FASTBIT_BIN\" -delimiter=\",\""

export INDEX_TYPE="FASTBIT:INDEX"
export INDEX_NAME_BASE="fastbit_idx"

RDI $INDEX_TYPE ${INDEX_NAME_BASE}_e.interval-equality_b.precision=2 $RDI_CMD -option=[e:interval-equality,b:precision=2]
