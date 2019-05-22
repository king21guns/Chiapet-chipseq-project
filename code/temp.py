    total = len(RAD)
    count = 0
    for chr in chrom_list:
        ups_value = []
        ins_value = []
        downs_value = []
        temp_RAD = RAD[RAD["chrom"] == chr]
        temp_file = file[file["chrom"] == chr]

        for i in xrange(len(temp_RAD)):
            if (count % 100 == 0) and (i != 0):
                print 'Dealing with %s: %d of %d\r' %(j,count,total),
                sys.stdout.flush()
            index = temp_RAD.index[i]
            start = RAD['start'][index]
            end = RAD['end'][index]
            a= temp_file[(temp_file['real_start']>=start-stream_num)&(temp_file['real_start']<=start)]
            if(len(a) != 0):
                ups_value.append(len(a))
            else:
                ups_value.append(0)
            b= temp_file[(temp_file['real_start']>=start)&(temp_file['real_start']<=end)]
            if(len(a) != 0):
                ins_value.append(len(b))
            else:
                ins_value.append(0)
            c= temp_file[(temp_file['real_start']>=end)&(temp_file['real_start']<=end+stream_num)]
            if(len(a) != 0):
                downs_value.append(len(c))
            else:
                downs_value.append(0)
            count += 1
        RAD['ups_'+histone][temp_RAD.index] = ups_value
        RAD['ins_'+histone][temp_RAD.index] = ins_value
        RAD['downs_'+histone][temp_RAD.index] = downs_value