"""
Created by Alejandro Plaza Cayon, 2020
"""

#Retrieve necessary modules
import requests, sys, os, time
from xml.etree import ElementTree as ElementTree

#Define a function posteriorly used to retrieve data references
def info_references(dato):
    try:
        if dato == 'disease':
            dat = com['description']
        elif dato == 'various':
            dat = com['text'][0]
        elif dato == 'catalytic':
            dat = com['reaction']
        if len(dat['evidences']) == 1 and dat['evidences'][0]['source']['url']:
            print(f"- REFERENCE: ", file=finf)
            print("> " + dat['evidences'][0]['source']['url'], file=finf)
        elif len(dat['evidences']) > 1:
            print(f"- REFERENCES: ", file=finf)
            for ref in dat['evidences']:
                try:
                    print("> " + ref['source']['url'], file=finf)
                except KeyError:
                    next
    except KeyError:
        next

#Function used for the retrieval of LSS variants, to make the process more flexible  
def try_write(data):
    if data != None:
        variants.write(str(data) + "\t")
    else:
        variants.write("\t")
        
#Dictionaries and variables needed for the code to work:
names = {"info":"txt", "curated_variants":"tsv"}       
ft = ['DOMAINS_AND_SITES', 'PTM', 'MUTAGENESIS', 'VARIANTS']
selected_species = {}
species_acc = {}
clinvar_ids = []
clinvar_acc = []
temp_ids = []
retstart = 0
total_count = 1

#Gene input:
gene_input = input("Enter a gene name: ")
gene_input = gene_input.upper().strip()

#Ask for desired species to query (human/mouse/rat)
sp = input(f"Look for human (Homo sapiens)? [y]/[n] ")
if sp == 'y':
    selected_species.update({'human': '9606'})
sp = input(f"Look for mouse (Mus musculus)? [y]/[n] ")
if sp == 'y':
    selected_species.update({'mouse': '10090'})
sp = input(f"Look for rat (Rattus norvegicus)? [y]/[n] ")
if sp == 'y':
    selected_species.update({'rat': '10116'})
print(f"\n>>> LOOKING FOR {gene_input}...\n")

#Now, the information in UniProt will be retrieved
for species, taxa in selected_species.items():
    #Accessing the database
    strc = {}
    td = {}
    requestURL = f"https://www.ebi.ac.uk/proteins/api/proteins?offset=0&size=100&reviewed=true&isoform=2&gene={gene_input}&taxid={taxa}"
    #Isoform= 0(just canonical), 1(isoforms), 2 (canonical and isoforms)
    #It searches by coding gene, it could search by protein but is more error-prone
    
    r = requests.get(requestURL, headers={ "Accept" : "application/json"})
    
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    
    #If there is no data, warn the user
    decoded = r.json()
    if not decoded:
        print(f"> No Uniprot information for {gene_input}")
    
    #Files creation and warnings for existing files
    else:
        for name, ext in names.items():
            if not os.path.isfile(f"./{gene_input}_{species}_{name}.{ext}"):
                if name == "info":
                    finf = open(rf"{gene_input}_{species}_info.txt","a")
                    finf.write(f"###UNIPROT INFORMATION FOR {gene_input}###\n\n")
                elif name == "curated_variants":
                    fres = open(rf"{gene_input}_{species}_curated_variants.tsv","a")
                    fres.write("Type\tPosition\tEffect\tReference\n")
            else:
                print(Warning(f">>> File ./{gene_input}_{species}_{name}.{ext} already exists!"))
                cont = input(">>> New information will be appended afterwards. Continue? [y]/[n]: ")
                cont = cont.lower()
                if cont == "y":
                    print("")
                    if name == "info":
                        finf = open(rf"{gene_input}_{species}_info.txt","a")
                        finf.write(f"\n\n###UNIPROT INFORMATION FOR {gene_input}###\n\n")
                    elif name == "curated_variants":
                        fres = open(rf"{gene_input}_{species}_curated_variants.tsv","a")
                        fres.write("Type\tPosition\tEffect\tReference\n")
                elif cont == "n":
                    print(">>> Aborting...")
                    try:
                        finf.close()
                        fres.close()
                        sys.exit()
                    except NameError:
                        sys.exit()
                else:
                    print(">>> Unknown command. Aborting...")
                    try:
                        finf.close()
                        fres.close()
                        sys.exit()
                    except NameError:
                        sys.exit()
                
    #Retrieving initial information: Name, accession  and entry IDs and organism:   
        for key, val in decoded[0].items():
            if key == 'accession':
                print("+ ACCESSION: " + val, file=finf)
                species_acc[species] = val
            elif key == 'secondaryAccession':
                if len(val) == 1:
                    print("+ SECONDARY ACCESSION: " + val[0], file=finf)
                elif len(val) > 1:
                    print("+ SECONDARY ACCESSIONS:", end = " ", file=finf)
                    for acc in val:
                        print(acc, end=" ", file=finf)
                    print("", file=finf)
            elif key =='id':
                print("+ UNIPROT ENTRY: " + val, file=finf)
                uniprot_entry = val
            elif key == 'organism':
                print("+ ORGANISM: " + val['names'][0]['value'] + " (Taxon " + 
                      str(val['taxonomy']) + ")", file=finf)
            elif key == 'protein':
                print("+ PROTEIN NAME: " + val['recommendedName']['fullName']['value'],
                      file=finf)
                
    #Retrieving the information associated to the comms list, as they follow the same data structure:
            elif key == 'comments':
                comms = ['FUNCTION', 'ACTIVITY_REGULATION', 'COFACTOR', 'SUBUNIT',
                         'SUBCELLULAR_LOCATION', 'TISSUE_SPECIFICITY', 'PATHWAY',
                         'DOMAIN', 'PTM', 'CAUTION', 'MISCELLANEOUS',
                         'SIMILARITY', 'POLYMORPHISMS']
                for com in val:
                    if com['type'] in comms:
                        data = com['type']
                        try:
                            print(f"\n+ {data}:", end = " ", file=finf)
                            text = com['text'][0]
                            print(text['value'], file=finf)
                            info_references('various')
                        except KeyError:
                            next
                        
    #Retrieve diseases information, which follows a different structure:
                    elif com['type'] == 'DISEASE':
                        print("\n+ DISEASE: ", end = " ", file=finf)
                        try:
                            print(com['diseaseId'].upper() + " - " + com['description']['value'],
                                  end = " ", file=finf)
                        except KeyError:
                            next
                            
                        try:
                            if len(com['text']) == 1:
                                print(f"\n- {gene_input} ROLE: ", end = " ", file=finf)
                                print(com['text'][0]['value'], file=finf)
                            elif len(com['text']) > 1:
                                print(f"\n- {gene_input} ROLE: ", end = " ", file=finf)
                                for role in com['text']:
                                    print(role['value'], file=finf)
                        except KeyError:
                            next
                        info_references('disease')
    
    #Retrieve cataliytic activity information:
                    elif com['type'] == 'CATALYTIC_ACTIVITY':
                        print("\n+ CATALYTIC ACTIVITY: ", end = " ", file=finf)
                        try:
                            print(com['reaction']['name'], file=finf)
                            info_references('catalytic')
                        except KeyError:
                            next
                            
    #Retrieve interactors:
                    elif com['type'] == 'INTERACTION':
                        print("\n+ INTERACIONS (Uniprot URLs): ", end = " ", file=finf)
                        for inter in com['interactions']:
                            print("\n> https://www.uniprot.org/uniprot/" + inter['accession2'], 
                                  end = " ", file=finf)
                            try:
                                print("- " + inter['gene'], end = " ", file=finf)
                            except KeyError:
                                next
                        print("\n", file=finf)
                            #Uncomment if you want to retrieve the name, organism and entry of the interactors:
                            #This could be updated to include further info about the interactors:
                            # new_requestURL = f"https://www.ebi.ac.uk/proteins/api/proteins/{inter['accession2']}"
                            # resp = requests.get(new_requestURL, headers={ "Accept" : "application/json"})
                            # if not resp.ok:
                            #     resp.raise_for_status()
                            #     sys.exit()
                            # new_decoded = resp.json()
                            # if not new_decoded:
                            #     print(f"No Uniprot information for {inter['accession2']}")
                            # else:
                            #     for k, v in new_decoded.items():
                            #         if k =='id':
                            #             print("-- UNIPROT ENTRY: " + v)
                            #         elif k == 'organism':
                            #             print("-- ORGANISM: " + v['names'][0]['value'] + " (Taxon " + 
                            #                   str(v['taxonomy']) + ")")
                            #         elif k == 'protein':
                            #             print("-- PROTEIN NAME: " + v['recommendedName']['fullName']['value'])
            
    #Retrieve sequence information:
            elif key == 'features':
                for each in val:
                    if each['type'] == 'CHAIN':
                        print("\n+ CHAIN LENGTH: " + each['end'] + " residues", file=finf)
                    elif each['category'] == 'TOPOLOGY':
                        if each['type'] == 'TRANSMEM':
                            td[f"{each['begin']}..{each['end']}"] = "Transmembrane - " + each['description']
                        if each['type'] == 'TOPO_DOM':
                            td[f"{each['begin']}..{each['end']}"] = each['description']
                        if each['type'] == 'INTRAMEM':
                            td[f"{each['begin']}..{each['end']}"] = "Intermembrane - " + each['description']
                    elif each['category'] == 'STRUCTURAL':
                        strc[f"{each['begin']}..{each['end']}"] = [each['type']]
                        if len(each['evidences']) == 1:
                            try:
                                strc[f"{each['begin']}..{each['end']}"].append(
                                    f"({each['evidences'][0]['source']['name']}" +
                                    f"- {each['evidences'][0]['source']['id']}" +
                                    f" - https://www.rcsb.org/structure/{each['evidences'][0]['source']['id']})")
                            except KeyError:
                                next
                        elif len(each['evidences']) > 1:
                            for ref in each['evidences']:
                                try:
                                    strc[f"{each['begin']}..{each['end']}"].append(
                                    f"({ref['source']['name']} - {ref['source']['id']}" +
                                    f"- https://www.rcsb.org/structure/{each['evidences'][0]['source']['id']})")
                                except KeyError:
                                    next
    #Write this retrieved information into the appropiate file:
                    if each['category'] in ft and each['type'] != 'VAR_SEQ':
                        try:                       
                            if each['type'] == 'VARIANT':
                                if each['begin'] == each['end']:
                                    fres.write(f"{each['type']}\t{each['alternativeSequence']}{each['begin']}\t{each['description']}\t")
                                else:
                                    fres.write(f"{each['type']}\t{each['alternativeSequence']}{each['begin']}-{each['end']}\t{each['description']}\t")
                            else:    
                                if each['begin'] == each['end']:
                                    fres.write(f"{each['type']}\t{each['begin']}\t{each['description']}\t")
                                else:
                                    fres.write(f"{each['type']}\t{each['begin']}-{each['end']}\t{each['description']}\t")
                            try:
                                if len(each['evidences']) == 1 and each['evidences'][0]['source']['url']:
                                    fres.write(f"{each['evidences'][0]['source']['url']}")
                                elif len(each['evidences']) > 1:
                                    for ref in each['evidences']:
                                        try:
                                            fres.write(f"{ref['source']['url']};")
                                        except KeyError:
                                            next
                            except KeyError:
                                next   
                            fres.write("\n")  
                        except KeyError:
                            next
        if td != {}:
            print("\n+ DOMAINS: ", file=finf)
            for k, v in td.items():
                print("> " + k + ": " + v, file=finf)
        if strc != {}:
            print("\n+ STRUCTURAL FEATURES: ", end = " ", file=finf)
            for k, v in strc.items():
                print("\n> " + k, end = " ", file=finf)
                for item in v:
                    print(item, end = " ", file=finf)
    
    
    #Retrieve the sequence from TOGOWS server, in fasta format, and write it in the same .txt file
    server = "http://togows.dbcls.jp/"
    ext = f"entry/uniprot/{uniprot_entry}.fasta"
     
    sequence = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
     
    if not r.ok:
      print(f"\n\n>>>UNABLE TO OBTAIN SEQUENCE FOR {gene_input}!!!!\n\n")
      sequence.raise_for_status()
      sys.exit()
    
    decoded_seq = sequence.text
    finf.write(f"\n\n+ {gene_input} CANNONICAL PROTEIN SEQUENCE, FASTA FORMAT:\n")
    print(f"{decoded_seq}", file=finf)    
    
    finf.close()
    fres.close()
    
print(">>> UniProt data retrieved, check for output files\n")

#With the initial selected species, variants from Large Scale Studies (LSS) will be retrieved:
for species, acc in species_acc.items():
    
    #Access the new server:
    requestURL = f"https://www.ebi.ac.uk/proteins/api/variation/{acc}"
    
    r = requests.get(requestURL, headers={ "Accept" : "application/json"})
    
    if not r.ok:
      r.raise_for_status()
      sys.exit()
    
    #Open a new file:
    if not os.path.isfile(f"./{gene_input}_{species}_LSS_variants.tsv"):
        variants = open(f"{gene_input}_{species}_LSS_variants.tsv","a")
    else:
        print(Warning(f">>> File ./{gene_input}_{species}_LSS_variants.tsv already exists!"))
        cont = input(">>> New information will be appended afterwards. Continue? [y]/[n]: ")
        cont = cont.lower()
        if cont == "y":
            print("")
            variants = open(f"{gene_input}_{species}_LSS_variants.tsv","a")
        elif cont == "n":
            print(">>> Aborting...")
            sys.exit()
        else:
            print(">>> Unknown command. Aborting...")
            sys.exit() 

    
    #Write the header of the file:
    variants.write("Begin\tEnd\tWT\tAlt\tCode\tMutation\tPolyPhen Score\tPolyPhen Pred.\t" + 
                   "SIFT Score\tSIFT Pred.\tID\tSource\tURL\tSource Type\n")
    responseBody = r.json()
    
    #Retrieve and write the data:
    for each in responseBody['features']:
        for ref in each['xrefs']:
            item = each.get('begin')
            try_write(item)
            item = each.get('end')
            try_write(item)
            item = each.get('wildType')
            try_write(item)
            item = each.get('alternativeSequence')
            try_write(item)
            if each['begin'] == each['end']:
                try_write(f"{each['wildType']}{each['begin']}{each['alternativeSequence']}")
            else:
                variants.write("\t")
            item = each.get('consequenceType')
            try_write(item)
            item = each.get('polyphenScore')
            try_write(item)
            item = each.get('polyphenPrediction')
            try_write(item)
            item = each.get('siftScore')
            try_write(item)
            item = each.get('siftPrediction')
            try_write(item)
            item = ref.get('id')
            try_write(item)
            item = ref.get('name')
            try_write(item)
            item = ref.get('url')
            if item != 'https://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?type=rs&rs=TCGA novel':
                try_write(item)
            elif item == 'https://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?type=rs&rs=TCGA novel':
                variants.write("\t")
            item = each.get('sourceType')
            try_write(item)
            variants.write("\n")      
    
    #Close the file                      
    variants.close()
print(">>> Large Scale Studies data retrieved, check for output file\n")

#Now ClinVar information will be retrieved.
#First, the variants' IDs will be identified for the selected gene.
#From the IDs, the accession codes will be retrieved, and finally, the final output info:
    
#The info ir retrieved in batches, to prevent the program from crashing:
while total_count > retstart:
    server = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    ext = f"esearch.fcgi?db=clinvar&term={gene_input}[gene]&retmax=500&retmode=json&retstart={retstart}"
     
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
     
    if not r.ok:
      r.raise_for_status()
      sys.exit()
     
    decoded = r.json()
    if decoded:
        total_count = int(decoded['esearchresult']['count'])
        retmax = int(decoded['esearchresult']['retmax'])
        temp_ids = decoded['esearchresult']['idlist']
        for ids in temp_ids:
            clinvar_ids.append(ids)
        retstart = retstart + retmax
    time.sleep(.34) #To prevent a ban from NCBI's API, there is a wait of 340ms between queries 
print(">>> ClinVar IDs retrieved\n")

for iden in clinvar_ids:
    ext = f"esummary.fcgi?db=clinvar&id={iden}&retmode=json"
     
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
     
    if not r.ok:
      r.raise_for_status()
      sys.exit()
     
    decoded = r.json()
    if decoded:
        try:
            clinvar_acc.append(decoded['result'][iden]['accession'])
        except Exception:
            next
    time.sleep(.34)
print(">>> ClinVar Accessions retrieved\n")

#Open a new file:
if not os.path.isfile(f"./ClinVar_{gene_input}.tsv"):
    fcv = open(f"ClinVar_{gene_input}.tsv", "a+")
else:
    print(Warning(f">>> File ./ClinVar_{gene_input}.tsv already exists!"))
    cont = input(">>> New information will be appended afterwards. Continue? [y]/[n]: ")
    cont = cont.lower()
    if cont == "y":
        print("")
        fcv = open(f"ClinVar_{gene_input}.tsv", "a+")
    elif cont == "n":
        print(">>> Aborting...")
        sys.exit()
    else:
        print(">>> Unknown command. Aborting...")
        sys.exit() 

#Write the file header
fcv.write("Clinvar ID\tVariant Type\tGene(s)\tHGVS\tAssembly\tPosition\tVariant Length\t" +
          "Genomic change\tEffect\tStudy\tSource type\tSource description\t" +
          "Source platform\tPhenotypes\tReference\n")
if clinvar_acc:
    for acc in clinvar_acc:
        ext = f"efetch.fcgi?db=clinvar&rettype=vcv&id={acc}&retmode=json"
         
        r = requests.get(server+ext, headers={ "Content-Type" : "application/xml"})
         
        if not r.ok:
          r.raise_for_status()
          sys.exit()
          
        methods = {}
        cites = {}
        phenotypes = []
        assoc_genes = []
        coords = {}
        var_len = {}
        alleles = {}
        study_op = None
        
        #Write the information from the XML file
        root = ElementTree.fromstring(r.content)
        if 'VariationID' in root[0].attrib:
            fcv.write(f"{root[0].attrib['VariationID']}\t")
        else:
            fcv.write("\t")
        if 'VariationType' in root[0].attrib:
            fcv.write(f"{root[0].attrib['VariationType']}\t")
        else:
            fcv.write("\t")
        for allele in root[0][2].findall("SimpleAllele"):
            for genelist in allele.findall("GeneList"):
                for gene in genelist.findall("Gene"):
                    assoc_genes.append(gene.attrib['Symbol'])
                fcv.write(f"{assoc_genes}")
                fcv.write("\t")
            for hgvs in allele.findall("Name"):
                fcv.write(f"{hgvs.text}\t")
            for loc in allele.findall("Location"):
                for x in loc:
                    if x.tag == 'SequenceLocation':
                        if 'innerStart' in x.attrib:
                            pos = x.attrib['innerStart'] + "-" + x.attrib['innerStop']
                            coords[x.attrib['Assembly']] = pos
                        elif 'start' in x.attrib:
                            pos = x.attrib['start'] + "-" + x.attrib['stop']
                            coords[x.attrib['Assembly']] = pos
                        if 'variantLength' in x.attrib:
                            var_len[x.attrib['Assembly']] = x.attrib['variantLength']
                        if 'referenceAlleleVCF' in x.attrib:
                            alleles[x.attrib['referenceAlleleVCF']] = x.attrib['alternateAlleleVCF']
                if 'GRCh38' in coords.keys():
                    fcv.write(f"GRCh38\t{coords['GRCh38']}\t")
                elif 'GRCh38' not in coords.keys() and 'GRCh37' in coords.keys():
                    fcv.write(f"GRCh37\t{coords['GRCh37']}\t")
                elif 'GRCh38' not in coords.keys() and 'GRCh37' not in coords.keys():
                    if len(coords) == 1:
                        fcv.write(f"{coords.keys()}\t{coord.values()}\t")
                    elif len(coords) > 1:
                        for k, v in coords.items():
                            fcv.write(f"{k}\t{v}\t")
                            break
                    elif len(coords) == 0:
                        fcv.write("\t\t")
                if 'GRCh38' in var_len.keys():
                    fcv.write(f"{var_len['GRCh38']}\t")
                elif 'GRCh38' not in var_len.keys() and 'GRCh37' in var_len.keys():
                    fcv.write(f"{var_len['GRCh37']}\t")
                elif 'GRCh38' not in var_len.keys() and 'GRCh37' not in var_len.keys():
                    if len(var_len) != 0:
                        for v in var_len.values():
                            fcv.write(f"{v}\t")
                    else:
                        fcv.write("\t")
                if len(alleles) != 0:
                    for k, v in alleles.items():
                        fcv.write(f"{k}>{v}\t")
                        break
                else:
                    fcv.write("\t")
        for clinical_l in root[0][2].findall('ClinicalAssertionList'):
            for clinical in clinical_l.findall('ClinicalAssertion'):
                for interp in clinical.findall('Interpretation'):
                    for descr in interp.findall('Description'):
                        if descr.text:
                            fcv.write(f"{descr.text}\t")
                        else:
                            fcv.write("\t")
                for study in clinical.findall('StudyDescription'):
                    study_op = study.text
                if study_op:
                    fcv.write(f"{study_op}\t")
                else:
                    fcv.write("\t")
                for citation in clinical.findall('Citation'):
                    for ids in citation.findall('ID'):
                        cites[ids.text] = ids.attrib['Source']
                    #print(citation[0].tag, citation[0].attrib)
                for obsev_l in clinical.findall('ObservedInList'):
                    for obsev in obsev_l.findall('ObservedIn'):
                        for method in obsev.findall('Method'):
                            for meth in method:
                                if meth.tag == 'Description':
                                    methods['description'] = meth.text
                                elif meth.tag == 'MethodType':
                                    methods['type'] = meth.text
                                elif meth.tag == 'TypePlatform':
                                    methods['platform'] = meth.text
                                if meth.tag == 'MethodAttribute':
                                    try:
                                        methods['attrib'] = meth[0].text
                                    except Exception:
                                        print("")
        for trait_l in root[0][2].findall('TraitMappingList'):
            for trait in trait_l.findall('TraitMapping'):
                if 'MappingValue' in trait.attrib:
                    if str(trait.attrib['MappingValue']).lower() != 'see cases':
                            phenotypes.append(trait.attrib['MappingValue'])
    
        if 'type' in methods.keys():
            fcv.write(f"{methods['type']}\t")
        else:
            fcv.write("\t")
            
        if 'description' in methods.keys() and 'attrib' in methods.keys():
            fcv.write(f"{methods['description']} - {methods['attrib']}\t")
        elif 'description' not in methods.keys() and 'attrib' in methods.keys():
            fcv.write(f"{methods['attrib']}\t")
        elif 'description' in methods.keys() and 'attrib' not in methods.keys():
            fcv.write(f"{methods['description']}\t")
        else:
            fcv.write("\t")
            
        if 'platform' in methods.keys():
            fcv.write(f"{methods['platform']}\t")
        else:
            fcv.write("\t")
            
        if len(phenotypes) != 0:
            if phenotypes[0] != 'not provided':
                fcv.write(f"{phenotypes}\t")
            else:
                fcv.write("\t")
            
        if len(cites) != 0:
            for i, j in cites.items():
                fcv.write(f"{j}: {i} ")
                fcv.write("\n")
        else:
            fcv.write("\n")
        time.sleep(.34)
        
    #Close ClinVar file:
    fcv.close()
    print(">>> ClinVar data retrieved, check for output file\n")
else:
    print(">>> There was a problem retrieving the accessions")
    print("The gene may not have associated data, or there may be a problem with the server")