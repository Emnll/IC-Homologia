/*A versão atual dos programas não possuem total tratamento de exceções pois se visou a estabilidade e confiabilidade acerca dos dados do mesmo, sendo possível, futuros trabalhos em cima da estrutura e confiabilidade do software.*/
Instruções para execução do programa:
	Os arquivos "kog.txt" e "kyva=gb.txt" devem estar na raiz do programa, com os programas : "cria_batch.py" e "busca_sequencia.py".

Executar na ordem : 
	1 - cria_batch;
	2 - busca_sequencia.

Variáveis Entrez:
	Modo anônimo: 
		Não informar nenhum email ou api key, simplesmente pressionando a tecla "Enter" para prosseguir.
	Modo autenticado: 
		Informar um e-mail válido e/ou a api key gerada pelo NCBI. 

Caso não tenha um email e uma api key: 

	Criar email : 
		Link abaixo para criar uma conta:
			account.ncbi.nlm.nih.gov/signup/?back_url=https%3A%2F%2Fwww.ncbi.nlm.nih.gov%2F
	API key: 
		Estando autenticado em sua conta:
			account settings > API Key Management 

Descrição do menu:
	cria_batch:
		1- cria diretório onde ficará os arquivos com a relação entre gi antigo e refseq atual.
		2- lê os arquivos do kog e pede a seleção do organismo previamente selecionado(hardcoded), realizar esta etapa múltiplas vezes(para cada organismo), afim de criar pastas para cada organismo contendo batches divididos em 199 entradas.

	*É obrigatório a primeira execução do 1º passo do menu
	**Ao menos uma execução do 2º passo se faz necessária para executar o segundo programa "busca_sequencia"
	***A 2ª etapa pode demorar!

	busca_sequencia:
		1- Criação diretório de sequências.
		2- Leitura dos batches para futura busca das sequências. Realizar este passo sempre antes de busca as sequências para definir qual pasta de organismo e qual batch será buscado.
		3- Cria diretório do organismo e escreve as sequências dos batches daquele organismo pretendido.
		4- Cria o diretório final para os resultados e junta todos os batches de sequência de cada organismo em um arquivo final de cada organismo.

	*É obrigatório a primeira execução do 1º passo do menu
	**As etapas dois e três ocorrem em conjunto, isso é, pode ser feito um ciclo das etapas para a confecção dos arquivos necessários demarcando os limites inferiores para futuros processamentos, assim evitando a execução única/intermitente do programa.
	***A 3ª etapa pode demorar!
	****Para utilizar as sequências em outros programas, como o orthofinder, é fortemente recomendado a realização da 4ª etapa para confecção de um único arquivo de cada organismo.
