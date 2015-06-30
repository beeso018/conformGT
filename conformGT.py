from optparse import OptionParser

parser = OptionParser()
parser.add_option("--ref", type=str)
parser.add_option("--test", type=str)
parser.add_option("--out", type=str)
options,args = parser.parse_args()

ref = open(options.ref, 'r')
test = open(options.test, 'r')
out = open(options.out+'.vcf', 'w')

class SNP(object):
    def __init__(self, ID, chrom, pos, ref, alt):
        self.id = ID
        self.chromosome = chrom
        self.position = pos
        self.reference = ref
        self.alternate = alt
        self.probs = []
    def addSample(self, sd, samp):
        self.probs.append(sd[samp])
    def getAnno(self):
        string = str(self.chromosome)+'\t'+ str(self.position)+'\t'+self.id+'\t'+self.reference+'\t'+self.alternate+'\t.\t.\t.\tGT:DS:GP'
        return string
    def getSNP(self):
        string = 'chr'+str(self.chromosome)+'\t'+ str(self.position)+'\t'+self.id+'\t'+self.reference+'\t'+self.alternate+'\t.\t.\t.\t.\tGT:DS:GP\t'+str(self.probs).replace('[','').replace(']','').replace("'","").replace(' ','\t').replace(',\t','\t')
        return string

def sampleDict(vcf):
    samples = {}
    for line in vcf:
        if line.startswith('#'):
            if line.startswith('#CHROM'):
                header = line.strip().split()
                for i in range(9, len(header)):
                    samples[i-9] = header[i]
                return samples

def snpDict(line, sampledict, gprobs=False):
    snps = {}
    ref = line[3]
    alt = line[4]
    for i in range(9, len(line)):
        info = line[i].strip().split(':')
        snp = info[0]
        if gprobs:
            gps = info[2].split(',')
        if snp == '0/0' or snp == '0|0':
            gt = ref+'/'+ref
            if gprobs:
                gp = gps[0]
        elif snp == '1/0' or snp == '0/1' or snp == '1|0' or snp == '0|1':
            gt = 'het'
            if gprobs:
                gp = gps[1]
        elif snp == '1/1' or snp == '1|1':
            gt = alt+'/'+alt
            if gprobs:
                gp = gps[2]
        elif snp == './.':
            gt = './.'
            gp = 0
        else:
            gt = './.'
            gp = 0
        if gprobs:
            snps[sampledict[i-9]] = str(gt)+':'+str(gp)
        else:
            snps[sampledict[i-9]] = gt
    return snps

def updateSnpDict(snpdict, o_r, o_a, n_r, n_a):
    for sample in snpdict:
        if snpdict[sample] == o_r+'/'+o_r:
            snpdict[sample] = n_a+'/'+n_a
        elif snpdict[sample] == 'het':
            snpdict[sample] = '0/1'
        elif snpdict[sample] == o_a+'/'+o_a:
            snpdict[sample] = n_r+'/'+n_r

def compareSnps(ref, test, o_r, o_a, n_r, n_a):
    na = 0
    het = 0
    href = n_r+'/'+n_r
    halt = n_a+'/'+n_a
    flip = False
    found = 0
    #    print(o_r, o_a, n_r, n_a)
    for key in test:
        try:
            tgt = test[key].split(':')
            if tgt[0] == href and ref[key] == halt:
                flip = True
                break
            elif tgt[0] == halt and ref[key] == href:
                flip = True
                break
            elif tgt[0] == ref[key]:
                if tgt[0] == 'het':
                    het += 1
                    found = 0
                flip = False
            elif test[key][:3] == './.' or ref[key] == './.':
                found = 0
                flip = False
            elif tgt[0] != ref[key]:
                if tgt[0] == o_r+'/'+o_r and ref[key] == href:
                    flip = False
                elif tgt[0] == o_a+'/'+o_a and ref[key] == halt:
                    flip = False
                else:
                    flip = True
                    break
            else:
                continue
        except KeyError:
            na += 1

    return flip

ref_samples = sampleDict(ref)
test_samples = sampleDict(test)

out.write('##fileformat=VCFv4.1\n##FORMAT=<ID=DS,Number=1,Type=Float,Description="estimated ALT dose [P(RA) + P(AA)]">\n##FORMAT=<ID=GP,Number=G,Type=Float,Description="Estimated Genotype Probability">\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
vcf_header = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'
for i in range(len(test_samples)):
    vcf_header += '\t'+test_samples[i]
out.write(vcf_header+'\n')

counter = 0
for line in test:
        counter += 1
        if counter%1000 == 0:
            print('Processed', counter, 'SNPs...')
        if line.startswith('#'):
            continue
        else:
            tl = line.strip().split()
            old_alt = tl[4]
            old_ref = tl[3]
            rl = [tl[0],0]
            try:
                found = False
                while rl[0] == tl[0] and int(rl[1]) <= int(tl[1]):
                        rl = ref.readline().strip().split()
                        if not rl == []:
                            if int(tl[1]) == int(rl[1]):
                                found = True
                                break
                if not found:
                    tl = line.strip().split()
                    if int(tl[1]) == int(rl[1]):
                        found = True
            except IndexError:
                continue
            if not found:
                continue
            elif found:
                new_ref = rl[3]
                new_alt = rl[4]
                ref_dict = snpDict(rl, ref_samples)
                test_dict = snpDict(tl, test_samples)
                s = SNP(rl[2], rl[0], rl[1], rl[3], rl[4])
                if compareSnps(ref_dict, test_dict, old_ref, old_alt, new_ref, new_alt):
                    updateSnpDict(test_dict, old_ref, old_alt, '0', '1')
                else:
                    updateSnpDict(test_dict, old_ref, old_alt, '1', '0')
                
                snps = ''
                for i in range(len(test_samples)):
                    snps += '\t'+test_dict[test_samples[i]]
                out.write(s.getAnno()+snps+'\n')
