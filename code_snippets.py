###Snippet 0###
#####################################

#import dependencies
import urllib2
import numpy as np
import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as sts

#import cobra


#####################################
###End###



###Snippet 1###
#####################################

#import bac_database table and print the names for the Fusobacterium strains

response = urllib2.urlopen('https://raw.githubusercontent.com/danielriosgarza/bacterial_passengers.py/master/input_files/bac_database.tsv')
bac_db = response.read()

for entry in bac_db.split('\n'):
  if 'Fusobacterium' in entry:
      a= entry.split('\t')
      print 'db_id: ', a[0], ' strain: ', a[1], ' NCBI_id: ', a[2],'\n'

    
#If you would like to download any figure or file generated with this notebook, first generate the file with 
#the data that you like to save and use something like the code below:

#from google.colab import files

#with open('example.txt', 'w') as f:
#  f.write('some content')

#files.download('example.txt')

#####################################
###End###



###Snippet 2###
#####################################

#get the Basal environment and store as a python_dict called mambo_env 
#containing [metabolite]=[uptake_rate]

response = urllib2.urlopen('https://github.com/danielriosgarza/bacterial_passengers.py/raw/master/results/mambo_env2.tsv')

md = response.read()

mambo_env={'EX_'+i.split('\t')[0]+'_e0':float(i.split('\t')[2]) for i in md.split('\n')[1:-1] }

#store metabolite names

code_to_name_met={'EX_'+i.split('\t')[0]+'_e0':i.split('\t')[1] for i in md.split('\n')[1:-1] }

print 'Basal environment dictionary\t', mambo_env

#####################################
###End###



###Snippet 3###
#####################################

#code below depends on cobrapy (critical dependencies are commented out)


def add_media(model, media_dict):
  '''  
  Add uptake constraints to exchange reactions
  
  input: 
  
  - python_dict [metbolite/transporter_id]:[uptake_limit]
  - cobra model
  
  outputs:
  
  -cobra model constrained by the uptake_limit
  -solution to the model
  '''
  
  for i in media_dict:
    try:#only take metabolites that the models can import 
      model.reactions.get_by_id(i).lower_bound = -abs(media_dict[i])
    except KeyError:
      pass
  
  for z in model.reactions:#set the bounds for all metabolite-uptakes that are not in the media_dict to a constant value of 1
    if 'EX_' in z.id and z.id not in media_dict:
      model.reactions.get_by_id(z.id).lower_bound = -1.
  
  solution = model.optimize()
  return solution.f

#with the dependecies installed, execute the code below:
#folder = '/models/'
#model_files= os.listdir(folder)

#mambo_obj={model.replace('.xml',''):add_media(cobra.io.read_sbml_model(folder+model), wes_diet) for model in model_files}

#####################################
###End###



###Snippet 4###
#####################################

#alternatively, take the pre computed solutions for the objective fluxes under the mambo environment

response = urllib2.urlopen('https://github.com/danielriosgarza/bacterial_passengers.py/raw/master/results/mambo_obj.tsv')

mambo_obj_f = response.read()

mambo_obj={i.split('\t')[0]:float(i.split('\t')[1]) for i in mambo_obj_f.split('\n')[1:-1]}


print 'mambo_obj\t', mambo_obj

#####################################
###End###



###Snippet 5###
#####################################

#Code used to construct table (partially comemented out because of dependencies). Takes several hours to complete on a regular pc or laptop.


def important_metabolites(model):
  '''
  obtain a list of important reactions
  
  '''
  s=model.optimize()
  sol=s.f
  important=[]
  for i in model.reactions:
    if '_e0' in i.id:
      cm = model.copy()
      cm.reactions.get_by_id(i.id).lower_bound=0
      s2 = cm.optimize()
      if s2.f<0.3*sol:
        important.append(i.id)
  return important

#with the dependecies installed, execute the code below:

#important_mets={}

#for model_name in model_files:
#    model = cobra.io.read_sbml_model(folder+model_name)
#    add_media(model,wes_diet)
#    important = important_metabolites(model)
#    important_mets[model_name.replace('xml','')]=np.array([1. if metabolite in important else 0. for metabolite in wes_diet])

#####################################
###End###



###Snippet 6###
#####################################

#obtain the precomputed data

response = urllib2.urlopen('https://github.com/danielriosgarza/bacterial_passengers.py/raw/master/results/IMPORTANT_metabolites_MAMBO_diet2.tsv')

important_f = response.read()

important_lines = important_f.split('\n')[1:-1]

metabolites=important_f.split('\n')[0].split('\t')[1:-1]

important_mets={}

for line in important_lines:
  a = line.split('\t')
  binary_vec = np.array(map(float, a[1:-1]))
  important_mets[a[0]]=binary_vec

#we must also take into account that diet metabolites could be absent from the model, so we must take a list of metabolites that are present in the models
#this will be used in the similarity score

response = urllib2.urlopen('https://github.com/danielriosgarza/bacterial_passengers.py/raw/master/input_files/metabolites_present.tsv')

present_f = response.read()

present_lines = present_f.split('\n')[1:-1]

metabolites_present=present_f.split('\n')[0].split('\t')[1:-1]

#obtain an index to convert metabolites_present to that of metabolites (above).
index_converter=np.array([metabolites.index(i) for i in metabolites_present])

#obtain a boolean vector for the presence or absence of metabolite in each model
present_mets={}

for line in present_lines:
  a = line.split('\t')
  boolean_vec = np.array(map(int, a[1:-1]), dtype=bool)[index_converter]#put them in the same order as in important_mets
  present_mets[a[0]]=boolean_vec
  

#show a simple plot with the distribution of important metabolites

important_met_dist=np.sum(important_mets.values(),axis=0)

met_sorter = np.argsort(important_met_dist)

plt.plot(important_met_dist[met_sorter], '-o'); plt.xlabel('metabolites'); plt.ylabel('bacteria') 

for i in xrange(len(metabolites)):
  print i, ' ', code_to_name_met[np.array(metabolites)[met_sorter][i]], ' ', important_met_dist[met_sorter][i]

#####################################
###End###



###Snippet 7###
#####################################

#Get the tumor-enriched metabolites
response = urllib2.urlopen('https://github.com/danielriosgarza/bacterial_passengers.py/raw/master/results/cancer_CM2.tsv')
cancer_C_m_f = response.read()

C_m={'EX_'+i.split('\t')[0]+'_e0':float(i.split('\t')[1]) for i in cancer_C_m_f.split('\n')[1:-1]}

print 'c_m ', C_m,'\n'

print 'Metabolites enriched in CRC tumor environment (according to 3 or more studies)\n'


for i in C_m:
  if C_m[i]==1.0:
    print code_to_name_met[i]

#index of cancer metabolites in the mambo_env metabolite list
can_met_index=np.array([metabolites.index(i) for i in C_m if C_m[i]==1.0])

#####################################
###End###



###Snippet 8###
#####################################

print '\tc\t(1)\t(0)\n\ni(n)\t(1)\t a\t b\n\n\t(0)\t c\t d'

#####################################
###End###



###Snippet 9###
#####################################

def binary_contigency(u, v):
  '''
  compute the entries of the contigency table for binary vectors (of the same size) u and v

  '''
  
  a=0
  b=0
  c=0
  d=0
  
  for i in xrange(len(u)):
    if u[i]==1 and v[i]==1:
      a+=1
    elif u[i]==1 and v[i]==0:
      b+=1
    elif u[i]==0 and v[i]==1:
      c+=1
    elif u[i]==0 and v[i]==0:
      d+=1
  
  return float(a),float(b),float(c),float(d)


def ochiai(u,v):
  '''
  Compute the Ochiai measure for binary vectors u and v
  
  '''
  
  a,b,c,d=binary_contigency(u,v)
  if a==0:
    return 0.
  else:
    return np.sqrt((a/(a+b))*(a/(a+c)))


#alternative binary measures of similarity (or distance may be used)

#def russellrao(u,v):
#  a,b,c,d=binary_contigency(u,v)
#  return a/(a+b+c+d)

S_n={}

#the index is consistent through the metabolites list

C_m_vec = np.array([C_m[i] for i in metabolites]) 

#present mets is used here
for bac in important_mets:
  S_n[bac]=ochiai(C_m_vec[present_mets[bac]],important_mets[bac][present_mets[bac]])

print 's(n) ', S_n

#to obtain a file with the computed s_n scores, uncoment the code below:

#f=file('s_scores.tsv', 'w')

#for strain_ in S_n:
#  f.write(strain_+'\t%.20f\n' %S_n[strain_])
#f.close()

#files.download('s_scores.tsv')

#####################################
###End###



###Snippet 10###
#####################################

#obtain lists with all strains and genera in the study
strains= np.array(S_n.keys())
genera=[]

for i in strains:
  if 'Candidatus' in i:
    name=i.split('_')
    genera.append(name[0]+'_'+name[1])
  else:
    genera.append(i.split('_')[0])

genera=np.array(genera)

#obtain list of genera enriched in CRC from github

response = urllib2.urlopen('https://raw.githubusercontent.com/danielriosgarza/bacterial_passengers.py/master/input_files/CRC_enriched_genera.tsv')

en_genera_f = response.read()

enriched_genera = en_genera_f.split('\n')[0:-1]


#rank the list of enriched genera by the average s_n

genera_len={}
others=[]

for i in xrange(len(genera)):
  if genera[i] in enriched_genera:
    if genera[i] not in genera_len:
      genera_len[genera[i]]=[S_n[strains[i]]]
    else:
      genera_len[genera[i]].append(S_n[strains[i]])
  else:
    others.append(S_n[strains[i]])
    
#get the mean s for enriched genera and others
mea={}

for i in genera_len:
  mea[i]=np.mean(genera_len[i])

others_mean = np.mean(others)  

#sort genera_by_mean
ek = np.array(mea.keys())
ev = np.array(mea.values())

sor = np.argsort(ev)
sor = sor[::-1]


#Define a color for each genera:


color_dict = {'Porphyromonas':'#fb9a99',

 'Prevotella':'red',
 'Gemella':'orange',
 'Fusobacterium': 'yellow',
 'Solobacterium':'cyan',
 'Anaerococcus':'blue',
 'Parvimonas':'lime',
 'Peptostreptococcus':'green',
 'Dialister':'#a6cee3',
 'Hungatella':'#33a02c',
 'Anaerotruncus':'#1f78b4',
 'Lachnoclostridium':'#cab2d6',
 'Blautia':'brown',
 'others':'#313233'}


#rank order (largest to smallest) based on the Sn vector
S_n_vec = np.array([S_n[i] for i in strains])

#Randomize possibly equal numbers 
S_n_diffs = np.array([i-j for i in S_n_vec for j in S_n_vec])#difference between values
min_pos_diff = min(S_n_diffs[S_n_diffs>0])#minimum positive difference

np.random.seed(666)

sorter = np.argsort(S_n_vec+np.random.uniform(0, min_pos_diff*0.9, len(S_n_vec)))
sorter=sorter[::-1]
rSn=S_n_vec[sorter]

#make the plots  

#colorplot
sorted_genera = genera[sorter]
sorted_strains = strains[sorter]

fig = plt.figure(0)
ax=fig.add_subplot(515)

for i in xrange(len(genera)):
  if sorted_genera[i] not in enriched_genera:#if from an enriched genera, plot a colored line
    ax.vlines(i, 0, 1, color=color_dict['others'],lw=.5)

for i in xrange(len(genera)):
  if sorted_genera[i] in enriched_genera:#if from an enriched genera, plot a colored line
    ax.vlines(i, 0, 1, color=color_dict[sorted_genera[i]],lw=1, label=genera[sorter][i])

plt.xlim(0,1561)

#boxplot

ax=fig.add_subplot(211)


#data for boxplots
data_boxplot=[genera_len[i] for i in ek[sor]]
data_boxplot.append(others)
data_boxplot=np.array(data_boxplot)

#data for labels
labels = [i for i in ek[sor]]
labels.append('non-CRC genera')
labels=np.array(labels)

#boxplot
medianprops = dict(linewidth=2.5, color='k')
#bplot = plt.boxplot(data_boxplot,notch=False, patch_artist=False, bootstrap=100000, labels=labels,medianprops=medianprops)

#scatter plot
for i in xrange(len(data_boxplot)-1):
  x = [i+1+np.random.uniform(-0.3, 0.3) for z in xrange(len(data_boxplot[i]))]
  plt.scatter(x, data_boxplot[i], c= color_dict[labels[i]],alpha=0.7)

  
x = [len(data_boxplot)+np.random.uniform(-0.3, 0.3) for z in xrange(len(data_boxplot[-1]))]

plt.scatter(x, data_boxplot[-1], c= color_dict['others'],alpha=0.7)

plt.xticks(rotation=90)

plt.ylabel('s-score')

plt.legend()


lgd = ax.legend(ncol=3 ,bbox_to_anchor=(1, 0.9),  prop={'size': 12})


#to dowload the figure uncomemment below:
#plt.savefig('Figure1.png',bbox_inches='tight', dpi=600)

#files.download('Figure1.png')

#####################################
###End###



###Snippet 11###
#####################################

def generate_D_null(n_permut, pos_set, complete_set):
  '''
  generate a null distribution for an accumulated weights 
  curve by randomizing the labels
  '''
  
  pos_set_=np.array(pos_set) 
  complete_set_=np.array(complete_set) 
  
  D_null=np.zeros(n_permut)
  Nall=float(len(complete_set))
  Npos = float(len([complete_set[i] for i in xrange(len(complete_set)) if complete_set[i] in pos_set]))
  
  L = (Npos-1)*np.sqrt((Nall-Npos)/Npos)
  for perm in xrange(n_permut):
    comp_set=complete_set_.copy()
    np.random.shuffle(comp_set)
    W=[0]
    for i in xrange(1,len(complete_set)):
      if comp_set[i] in pos_set:
        W.append(W[-1]+ np.sqrt((Nall-Npos)/Npos))
      else:
        W.append(W[-1]-np.sqrt(Npos/(Nall-Npos)))
    D_null[perm]=max(W)
  return np.array(D_null)/L


def generate_D_null_rpos(n_permut, size_pos, complete_set):
  '''
  generate a null distribution for an accumulated 
  weights curve by randomizing the positive set
  
  '''
  
  D_null=np.zeros(n_permut)
  Nall =float(len(complete_set))
  comp_s=list(set(complete_set))
  weights= np.array([complete_set.count(i) for i in comp_s])
  weights = weights/float(sum(weights))
  
  for perm in xrange(n_permut):
    pos_set = np.random.choice(comp_s, size_pos, replace=1,p=weights)
    Npos = float(len([complete_set[i] for i in xrange(len(complete_set)) if complete_set[i] in pos_set]))
    L = (Npos-1)*np.sqrt((Nall-Npos)/Npos)
    W=[0]
    for i in xrange(1,len(complete_set)):
      if complete_set[i] in pos_set:
        W.append(W[-1]+ np.sqrt((Nall-Npos)/Npos))
      else:
        W.append(W[-1]-np.sqrt(Npos/(Nall-Npos)))
    D_null[perm]=max(W)/L
  return np.array(D_null)




#Obtain Nall, Npos, and L

Nall=float(len(rSn))

Npos = float(len([rSn[i] for i in xrange(len(rSn)) if genera[sorter][i] in enriched_genera]))

L = (Npos-1)*np.sqrt((Nall-Npos)/Npos)

#Build the weight curve
W=[0]

for i in xrange(1,len(strains)):
  if genera[sorter][i] in enriched_genera:
    W.append(W[-1]+ np.sqrt((Nall-Npos)/Npos))
  else:
    W.append(W[-1]-np.sqrt(Npos/(Nall-Npos)))

#constrain to the maximum of 1              
W = np.array(W)/L             


D_null = generate_D_null(10000, enriched_genera, genera)
D_null_rpos = generate_D_null_rpos(10000, 13,list(sorted_genera))


emp_pv = (1.+len([i for i in D_null if i>=max(W)]))/(len(D_null))
emp_pv_rpos = (1.+len([i for i in D_null_rpos if i>=max(W)]))/(len(D_null_rpos))

z_score = (max(W)-np.mean(D_null))/np.std(D_null)
z_score_rpos = (max(W)-np.mean(D_null_rpos))/np.std(D_null_rpos)

pv = sts.norm.sf(abs(z_score))
pv_rpos = sts.norm.sf(abs(z_score_rpos))

#####################################
###End###



###Snippet 12###
#####################################

print 'L = %.7f' %L

print 'max(W) = %.7f' %max(W)

print 'mean(Dnull) = %.7f' %np.mean(D_null)

print 'mean(Dnull_rpos) = %.7f' %np.mean(D_null_rpos)

print 'empirical_pvalue = %.7f' %emp_pv

print 'empirical_pvalue_rpos = %.7f' %emp_pv_rpos

print 'estimated pvalue = %.7f' %pv

print 'estimated pvalue_rpos = %.7f' %pv_rpos
           

#plot the accumulated weight curve and the compariso of Dnull and max(W)  

f, axarr = plt.subplots(2,2)
axarr[0,1].hist(D_null, 50, density=1, alpha=0.5, color='r')
axarr[0,1].hist(D_null, 50, density=1, histtype='step', lw=1.5, color='r')
axarr[0,1].set_xticks([np.mean(D_null)])
axarr[0,1].set_xticklabels(['%.2f'%np.mean(D_null)])
axarr[0,1].vlines(max(W),0,15., lw=2,color='g')

axarr[1,1].hist(D_null_rpos, 50, density=1, alpha=0.5, color='orange')
axarr[1,1].hist(D_null_rpos, 50, density=1, histtype='step', lw=1.5, color='orange')
axarr[1,1].set_xticks([np.mean(D_null_rpos), max(W)])
axarr[1,1].set_xticklabels(['%.2f'%np.mean(D_null_rpos), '%.2f\nmax(W)'%max(W)])
axarr[1,1].vlines(max(W),0,20., lw=2,color='g')


  
plt.subplot(121).plot(W, color=color_dict['others'], lw=1)
for i in xrange(len(genera[sorter])):
    if genera[sorter][i] in enriched_genera:
        plt.subplot(121).scatter(i, W[i], c=color_dict[genera[sorter][i]])
    
        
plt.ylabel('cumulative weight distribution (W)')

#plt.savefig('Figure2.png',bbox_inches='tight',dpi=600)

#files.download('Figure2.png')

#####################################
###End###



###Snippet 13###
#####################################

#Compute fn_null

def null_d(model, media_dict, n_mets, n_points):
  '''
  generate (fn)null
  
  '''
  add_media(model, media_dict)
  s=model.optimize()
  sol=s.f
  k=np.array(media_dict.keys())
  rfc=np.zeros(n_points)
  for i in xrange(n_points):
    l = k.copy()
    mc = model.copy()
    np.random.shuffle(l)
    
    d={z:10000 for z in l[0:n_mets]}
    add_media(mc, d)
    s2 = mc.optimize()
    rfc[i]=np.log2(s2.f/sol)
    return rfc
  
#comented out because of dependency with cobrapy 

#for mod in folder:
#  model = cobra.io.read_sbml_model(folder+mod)
#  sampling_distribution =  null_d(model, media)
#  f = file(path_to_output+mod.replace('.xml', '.txt'),'w')
#  for i in sampling_distribution:
#    f.write('%.20f\n'%i)
#    f.close()
    
#####################################
###End###



###Snippet 14###
#####################################

#Obtaining bn from github
response = urllib2.urlopen('https://github.com/danielriosgarza/bacterial_passengers.py/raw/master/results/bn_scores_MAMBO.tsv')

Bn_f = response.read()

B_n = {i.split('\t')[0]:float(i.split('\t')[1]) for i in Bn_f.split('\n')[0:-1]}

print 'bn ', B_n

#####################################
###End###



###Snippet 15###
#####################################

#rank the list of enriched genera by the average B_n

genera_len_b={}
others_b=[]

for i in xrange(len(genera)):
  if genera[i] in enriched_genera:
    if genera[i] not in genera_len_b:
      genera_len_b[genera[i]]=[B_n[strains[i]]]
    else:
      genera_len_b[genera[i]].append(B_n[strains[i]])
  else:
    others_b.append(B_n[strains[i]])
    
#get the mean s for enriched genera and others
mea_b={}

for i in genera_len_b:
  mea_b[i]=np.mean(genera_len_b[i])

others_mean_b = np.mean(others_b)  


#rank order (largest to smallest) based on the Bn vector
B_n_vec = np.array([B_n[i] for i in strains])

#Randomize possibly equal numbers 
B_n_diffs = np.array([i-j for i in B_n_vec for j in B_n_vec])#difference between values
min_pos_diff = min(B_n_diffs[B_n_diffs>0])#minimum positive difference

np.random.seed(666)

sorter_b = np.argsort(B_n_vec+np.random.uniform(0, min_pos_diff*0.9, len(B_n_vec)))
sorter_b=sorter_b[::-1]
rBn=B_n_vec[sorter_b]

#make the plots  

#colorplot
sorted_genera_b = genera[sorter_b]
sorted_strains_b = strains[sorter_b]

fig = plt.figure(0)
ax=fig.add_subplot(515)

for i in xrange(len(genera)):
  if sorted_genera_b[i] not in enriched_genera:
    ax.vlines(i, 0, 1, color=color_dict['others'],lw=.5)

for i in xrange(len(genera)):
  if sorted_genera_b[i] in enriched_genera:
    ax.vlines(i, 0, 1, color=color_dict[sorted_genera_b[i]],lw=1, label=genera[sorter_b][i])

plt.xlim(0,1561)

#boxplot

ax=fig.add_subplot(211)

#data for boxplots
data_boxplot=[genera_len_b[i] for i in ek[sor]]
data_boxplot.append(others_b)
data_boxplot=np.array(data_boxplot)

#boxplot
medianprops = dict(linewidth=2.5, color='k')
#bplot = plt.boxplot(data_boxplot,notch=False, patch_artist=False, bootstrap=100000, labels=labels,medianprops=medianprops)

#scatter plot
for i in xrange(len(data_boxplot)-1):
  x = [i+1+np.random.uniform(-0.3, 0.3) for z in xrange(len(data_boxplot[i]))]
  plt.scatter(x, data_boxplot[i], c= color_dict[labels[i]],alpha=0.7)

  
x = [len(data_boxplot)+np.random.uniform(-0.3, 0.3) for z in xrange(len(data_boxplot[-1]))]

plt.scatter(x, data_boxplot[-1], c= color_dict['others'],alpha=0.7)

plt.xticks(rotation=90)

plt.ylabel('s-score')

plt.legend()


lgd = ax.legend(ncol=3 ,bbox_to_anchor=(1, 0.9),  prop={'size': 12})


#to dowload the figure uncomemment below:
#plt.savefig('Figure1.png',bbox_inches='tight', dpi=600)

#files.download('Figure1.png')

#####################################
###End###


###Snippet 16###
#####################################

#Obtain Nall, Npos, and L

Nall=float(len(rBn))

Npos = float(len([rBn[i] for i in xrange(len(rBn)) if genera[sorter_b][i] in enriched_genera]))

L = (Npos-1)*np.sqrt((Nall-Npos)/Npos)

#Build the weight curve
W=[0]

for i in xrange(1,len(strains)):
  if sorted_genera_b in enriched_genera:
    W.append(W[-1]+ np.sqrt((Nall-Npos)/Npos))
  else:
    W.append(W[-1]-np.sqrt(Npos/(Nall-Npos)))

#constrain to the maximum of 1              
W = np.array(W)/L             


D_null = generate_D_null(10000, enriched_genera, genera)
D_null_rpos = generate_D_null_rpos(10000, 13,list(sorted_genera_b))


emp_pv = (1.+len([i for i in D_null if i>=max(W)]))/(len(D_null))
emp_pv_rpos = (1.+len([i for i in D_null_rpos if i>=max(W)]))/(len(D_null_rpos))

z_score = (max(W)-np.mean(D_null))/np.std(D_null)
z_score_rpos = (max(W)-np.mean(D_null_rpos))/np.std(D_null_rpos)

pv = sts.norm.sf(abs(z_score))
pv_rpos = sts.norm.sf(abs(z_score_rpos))

#####################################
###End###
