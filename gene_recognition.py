############################
#   제출자 1771089 김소현  #
#   2019년 06월 22일       #
############################
#idle입력을 조금씩 복사 붙여넣기 한 것으로 누락된 내용이나
#옮기는 과정에서 indentation이 틀린 부분이 있을 수 있습니다.


import nltk
from nltk.corpus import *
import wikipedia
import random
import re

stopwordslist=stopwords.words('english')
def hasNumbers(inputString):
	return any(char.isdigit() for char in inputString)



################ 아래 def가 알수없는 indentation 에러가 발생합니다. ###################
#실행당시에는 아무 문제 없었고, 또 복사 붙여넣기를 해서 실행하면 아무런 에러가 발생하지
#않습니다. 유독 이 창에서 run module을 선택하면 에러가 발생(ㅠㅠ) 하니
#혹시 테스트 해 보실거라면 여기서부터 아래는 번거로우시겠지만
#복사 붙여넣기를 해 주세요.

def feature_extracter(word):
	features={}


	#pos tag결과
	pos=nltk.pos_tag([word])
	if (pos==[(word,'NN')]):
		features['NN']=True
		features['DT']=False
	elif(pos==[(word,'DT')]):
	      features['NN']=False
	      features['DT']=True
	else:
	      features['NN']=False
	      features['DT']=False

	#그 밖
	features['contain number']=hasNumbers(word)
	features['isupper']=word.isupper()
	popular=['TP53', 'TNF', 'EGFR', 'VEGFA', 'APOE', 'IL6', 'TGFB1', 'MTHFR', 'ESR1','AKT1']
	features['contains popular words']=False
	for token in popular:
		if(not features['contains popular words']):
			features['contains popular words']= (token in word)
	features['stopwords']=word in stopwordslist

	noresult=0
	try:
		#위키피디아 사용
		ny=wikipedia.page(word)
		wikicontent=ny.content
		#features['wiki-proteni?']='protein' in wikicontent
		#features['wiki-genes?']='gene' in wikicontent

	except wikipedia.exceptions.DisambiguationError as e:

                try:
                        s=random.choice(e.options)
                        ny=wikipedia.page(s)
                        wikicontent=ny.content
                except wikipedia.exceptions.DisambiguationError as e:
                        features['wiki-proteni?']=False
                        features['wiki-genes?']=False
                        features['applied by gene?']=False
                        noresult=1 #goto문 대신 사용할 flag
                except wikipedia.exceptions.WikipediaException as f:
                        features['wiki-proteni?']=False
                        features['wiki-genes?']=False
                        features['applied by gene?']=False
                        noresult=1 #goto문 대신 사용할 flag

	except wikipedia.PageError as ee:
                features['wiki-proteni?']=False
                features['wiki-genes?']=False
                features['applied by gene?']=False
                noresult=1 #goto문 대신 사용할 flag

	except wikipedia.exceptions.WikipediaException as eee:

                features['wiki-proteni?']=False
                features['wiki-genes?']=False
                features['applied by gene?']=False
                noresult=1 #goto문 대신 사용할 flag

	if (noresult !=1):
		try:
			features['wiki-proteni?']='protein' in wikicontent
			features['wiki-genes?']='gene' in wikicontent
			text=nltk.Text(nltk.word_tokenize(wikicontent)) #래핑
			bigram_wiki=list(nltk.bigrams(text))
			word_bf_gene = [a[0] for a in bigram_wiki if a[1]=='gene']
			if(word in word_bf_gene):
				features['applied by gene?']=True
			else:
				features['applied by gene?']=False
		except wikipedia.exceptions as f:
			features['wiki-proteni?']=False
			features['wiki-genes?']=False
			features['applied by gene?']=False
	
	return features


#예시 feature extractoer
#feature_extracter('RTNF')


#training data 만들기############################################################
#positive data#
corpusroot="C:\\Users\\yorco\\Desktop\\imp\\" #파일의 경로(필요시 수정)
nc=PlaintextCorpusReader(corpusroot,'goa_human.gaf',encoding='utf-8')
raw=nc.raw()
genelist=re.findall('UniProtKB	([A-Za-z0-9]+)	[A-Za-z0-9-:]+',raw)
genelist=set(genelist)
genelist_t=re.findall('\t[CPF]{1}\t(.+?)\s[A-Z]+',raw)
genelist_t=set(genelist_t)
genelist=list(genelist)+list(genelist_t)
#len(genelist) #길이 확인, 약 3만개 일시 정상작동함
random.shuffle(genelist)
#genelist=genelist[:100] 데이터 양을 가볍게 줄이고 싶을 때
proteinlist=[ (feature_extracter(token),'YES') for token in genelist] #warning발생하지만 잘 작동됨


#negative data#
falselist=stopwords.words('english')
falselist.extend(treebank.words('wsj_0003.mrg'))
falselist=set(falselist)
falselist=list(falselist)
webstring=webtext.raw()
webt=nltk.sent_tokenize(webstring)
for i in range(20):
	falselist.append(webt[i+5])
random.shuffle(falselist)
#falselist=falselist[:100] 데이터의 양을 가볍게 줄이고 싶을 때
notgenes=[(feature_extracter(token),'No') for token in falselist]
featureset=proteinlist+notgenes
random.shuffle(featureset)
index=int(len(featureset)*0.9)
train_set=featureset[:index]
test_set=featureset[index:]
classifier=nltk.NaiveBayesClassifier.train(train_set)
print(nltk.classify.accuracy(classifier,test_set))
classifier.show_most_informative_features()

#test#
#classifier.classify(feature_extracter('P35'))
#classifier.classify(feature_extracter('hello-nice to meet you'))
#classifier.classify(feature_extracter('APPLE'))
#classifier.classify(feature_extracter('APPLE2'))

#gnitest############################################################################
croot="C:\\Users\\yorco\\Desktop\\imp\\gni" #경로 : 필요시 수정
cp=PlaintextCorpusReader(croot,'.*\.txt',encoding='utf-8')
gniwords=cp.words()
#gniwords=gniwords[:100]#테스트용 소량
gnigene=[]
for i in gniwords:
	result=classifier.classify(feature_extracter(i))
	if (result=='YES'):
		gnigene.append(i)
