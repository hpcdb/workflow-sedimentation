#Caminho até o jar do RDI
export RDI_DIR='/home/thaylon/Commands/RDI-1.0-jar-with-dependencies.jar'

#Criando o comando "RDI"
alias RDI="java -jar \"$RDI_DIR\""

#Caminho até diretório bin/ do PostgresRAW
export BIN='/home/thaylon/Softwares/PostgresRAW/bin/'

#Caminho até diretório pgData
export PG_DATA='/home/thaylon/pgData'

#As colunas presente no arquivo CSV
#Existem 3 tipos de colunas, sendo eles
#Arquivo -> file
#Texto -> text
#Número inteiro -> number
#Número decimal -> numeric:<número de casas depois da vírgula>
export COLUMNS='[u:numeric:50,v:numeric:50,w:numeric:50,p:numeric:50,s:numeric:50,d:numeric:50,vtkValidPointMask:numeric,arc_length:numeric:50,Points0:numeric:50,Points1:numeric:50,Points2:numeric:50]'

#Indexação por regex, no caso todos arquivos .csv
export RDI_CMD="POSTGRES_RAW:INDEX di_csv_idx_postgres_raw $PWD \"*.csv\" $COLUMNS -bin=\"$BIN\" -delimiter=\",\" -pgData=\"$PG_DATA\" -option=[header:true]"

RDI $RDI_CMD
