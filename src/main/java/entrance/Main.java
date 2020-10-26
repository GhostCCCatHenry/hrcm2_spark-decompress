package entrance;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileStatus;
import org.apache.hadoop.fs.FileSystem;
import org.apache.hadoop.fs.Path;
import util.MatchEntry;

import java.io.*;
import java.net.URI;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static java.lang.Math.ceil;

public class Main {

    private static final int MAX_SEQ_NUM = 2000;//maximum sequence number
    private static final int MAX_CHA_NUM = 1 << 28;//maximum length of a chromosome
    private static final int LINE_CHA_NUM =200;
    private static final int PERCENT = 15; //the percentage of compressed sequence number uses as reference
    private static final int kMerLen = 13; //the length of k-mer
    private static final int kmer_bit_num = 2 * kMerLen; //bit numbers of k-mer
    private static final int hashTableLen = 1 << kmer_bit_num; // length of hash table
    private static final int VEC_SIZE = 1 <<20; //length for other character arrays
    private static final int min_rep_len = 15;   //minimum replace length, matched string length exceeds min_rep_len, saved as matched information
    private static long startTime;

    private static String identifier;
    private static int lineWidth, ref_code_len, seq_code_len, ref_low_len, seq_low_len, diff_low_len, nCha_len, spe_cha_len, seqNumber, seqBucketLen;
    private static int sec_seq_num; //the referenced sequence number used for second compress
    private static final char []ref_code=new char[MAX_CHA_NUM];
    private static final char []seq_code=new char[MAX_CHA_NUM];//mismatched subsequence
    private static final int []low_loc= new int[VEC_SIZE/2]; //lowercase tuple location

    private static int[] ref_low_begin;
    private static int[] ref_low_length;
    private static int[] seq_low_begin;
    private static int[] seq_low_length;
    private static int[] diff_low_begin;
    private static int[] diff_low_length;
    private static int[] nCha_begin;
    private static int[] nCha_length;
    private static int[] spe_cha_pos;
    private static int[] spe_cha_ch;

    private static final List<String> seqPath = new ArrayList<>();
    private static final List<String> seqName = new ArrayList<>();
//    private static int[] lineWidth_vec = new int[2000];
    private static List <MatchEntry> matchResult;                 //存储一个序列的初次匹配结果
    private static List <MatchEntry> misMatchEntry;               //存储二次匹配时未匹配matchentry
    private static List <List <MatchEntry> > matchResult_vec;
    private static List <List<Integer> > seqLoc_vec;
    private static List<int []> seqBucket_vec;

    // allocate memory
    private static void initial() {
        ref_low_begin = new int[VEC_SIZE];
        ref_low_length = new int[VEC_SIZE];
        seq_low_begin = new int[VEC_SIZE];
        seq_low_length = new int[VEC_SIZE];
        diff_low_begin = new int[VEC_SIZE];
        diff_low_length = new int[VEC_SIZE];
        nCha_begin = new int[VEC_SIZE];
        nCha_length = new int[VEC_SIZE];
        spe_cha_ch = new int[VEC_SIZE/2];
        spe_cha_pos = new int[VEC_SIZE/2];

//        identifier_vec = new ArrayList<>();
        matchResult = new ArrayList<>();
        misMatchEntry = new ArrayList<>();
        matchResult_vec=new ArrayList<>();
        seqLoc_vec=new ArrayList<>();
        seqBucket_vec=new ArrayList<>();
        //lineWidth_vec.reverse();
        //seqBucket_vec.reserve(seqNumber);
        //seqLoc_vec.reserve(seqNumber);
    }

    private static final char []integerEncoding = { 'A', 'C', 'G', 'T' };

    private static void readFile(String filePath) throws NullPointerException{
        File fp =new File(filePath);
        String[] names = fp.list();
        if(names!=null)
            Arrays.sort(names);//这里添加一个排序操作，对names数组进行排序，然后再把它放进seqPath数据组里去
            for (String a: names) {
                if(!a.contains(".crc")&&a.contains("part-0"))
                    seqPath.add(a);
               // System.out.println(a);
            }
        seqNumber = seqPath.size();
    }
    private static void  readName(String namePath) {

        try {
            String name = null;
            Configuration conf = new Configuration();
            conf.set("fs.hdfs.impl", "org.apache.hadoop.hdfs.DistributedFileSystem");
            FileSystem fs = FileSystem.get(new URI("hdfs://master:9000"), conf);
            FileStatus[] fss = fs.listStatus(new Path(namePath));
            for(FileStatus f:fss){
                seqName.add(f.getPath().getName());
            }

        } catch (IOException | URISyntaxException e) {
            e.printStackTrace();
        }
        seqNumber = seqName.size();
    }

    private static void referenceSequenceExtraction(String str_referenceName) {
        int _seq_code_len = 0, _ref_low_len = 1, letters_len = 0;//record lowercase from 1, diff_lowercase_loc[i]=0 means mismatching
        char temp_cha;
        boolean flag = true;
        String str;
        char[] cha;      //the content of one line
        File fp=new File(str_referenceName);
        BufferedReader br;

        try {
            br=new BufferedReader(new InputStreamReader(new DataInputStream(new FileInputStream(fp))));
            br.readLine();
            while ((str=br.readLine())!=null)
            {
                cha=str.toCharArray();
                for (int i=0;i<cha.length;i++){
                    temp_cha=cha[i];
                    if(Character.isLowerCase(temp_cha)){
                        if (flag) //previous is upper case
                        {
                            flag = false; //change status of flag
                            ref_low_begin[_ref_low_len] = letters_len;
                            letters_len = 0;
                        }
                        temp_cha = Character.toUpperCase(temp_cha);
                    }
                    else {
                        if (!flag)  //previous is lower case
                        {
                            flag = true;
                            ref_low_length[_ref_low_len++] = letters_len;
                            letters_len = 0;
                        }
                    }
                    if (temp_cha == 'A' || temp_cha == 'C' || temp_cha == 'G' || temp_cha == 'T')
                        ref_code[_seq_code_len++] = temp_cha;
                    letters_len++;
                }
            }
            br.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        if (!flag)  //if flag=false, don't forget record the length
            ref_low_length[_ref_low_len++] = letters_len;

        ref_code_len = _seq_code_len;
        ref_low_len = _ref_low_len - 1;
//        System.out.println("Extraction of reference sequence complete. Reference code length: %d. Lowercase length: %d.\n"+ ref_code_len+"\t"+ref_low_len);

    }

    private static void seqLowercaseReading(int _seq_low_len) {
        int loc;
        //printf("The seq_low_len is %d, the diff_low_len is %d\n", _seq_low_len, _diff_low_len);
        for (int i = 0, j =0; i < _seq_low_len; i++) {
            loc = low_loc[i];
            if (loc == 0) {
                seq_low_begin[i] = diff_low_begin[j];
                seq_low_length[i] = diff_low_length[j++];
            }
            else {
                seq_low_begin[i] = ref_low_begin[loc];
                seq_low_length[i] = ref_low_length[loc];
            }
        }
    }

    private static int[] runLengthDecoding(BufferedReader br, int tolerance) {
//        File fp = new File(filename);
        int code_len, temp;
        int[] vec;
        List<Integer> code = new ArrayList<>();
        try {
//            BufferedReader br = new BufferedReader(new FileReader(fp));
            code_len = Integer.parseInt(br.readLine());
            for (int i = 0; i < code_len; i++)
            {
//                fscanf(fp, "%d", &temp);
                temp = Integer.parseInt(br.readLine());
                code.add(temp);
            }

            int length = 0;
            for (int i = 1; i < code_len; i += 2)
                length += code.get(i);
            //printf("The length of sequence lowercase length is : %d\n", length);
//            seq_low_len = length;
            if (length > 0) {
                vec = new int[length];
                int k = 0;
                for (int i = 0; i < code_len; i += 2)
                    for (int j = 0; j < code.get(i+1); j++)
                        vec[k++] = code.get(i)+j*tolerance;
                //printf("The length of sequence lowercase length is : %d\n", k);
                return vec;
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
//        fscanf(fp, "%d", &code_len);
        return null;
    }

    private static void readIdentifierData(BufferedReader br, List<String> vec) {
        int code_len, temp_int;
        String str;
        List<Integer> code = new ArrayList<>();
//        char []temp_str = new char[LINE_CHA_NUM];
        try {
            code_len = Integer.parseInt(br.readLine());
//            fscanf(fp, "%d", &code_len);
            for (int i = 0; i < code_len; i++) {
                temp_int = Integer.parseInt(br.readLine());
//                fscanf(fp, "%d", &temp_int);
                code.add(temp_int);
            }
            for (int i = 0; i < code_len; i++) {
                while ((str=br.readLine())!=null) {
                    if(str.charAt(0) == '\n') continue;
//                    temp_str = str.toCharArray();
//                    if (temp_str[0] == '\n') continue;
                    identifier = str;
                    for (int j = 0; j < code.get(i); j++)
                        vec.add(identifier);
                    break;
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void readPositionRangeData(BufferedReader br, int _vec_len, int[] _vec_begin, int[] _vec_length) {
//        fscanf(fp, "%d", &_vec_len);
        try {
            for (int i = 0; i < _vec_len; i++){
                _vec_begin[i] = Integer.parseInt(br.readLine());
                _vec_length[i] = Integer.parseInt(br.readLine());
            }

        } catch (IOException e) {
            e.printStackTrace();
        }

//            fscanf(fp, "%d%d", &_vec[i].begin, &_vec[i].length);
    }

/*    private static void readSpeChaData(FILE *fp, int &_vec_len, POSITION_SPE_CHA *&_vec)
    {
        _vec = new POSITION_SPE_CHA[_vec_len];
        for (int i = 0; i < _vec_len; i++)
            fscanf(fp, "%d", &_vec[i].pos);

        //read special character information
        vector<int> arr;
        int size, temp;
        fscanf(fp, "%d", &size);
        for (int i = 0; i < size; i++)
        {
            fscanf(fp, "%d", &temp);
            arr.push_back(temp);
        }
        if (size != 1)
        {
            int bit_num = ceil(log(size) / log(2));
            int v_num = floor(32.0 / bit_num);

            for (int i = 0; i < _vec_len; )
            {
                unsigned int v;
                fscanf(fp, "%u", &v);
                vector<int> temp_arr;

                int temp_i = i;
                for (int j = 0; j < v_num && temp_i < _vec_len; j++, temp_i++)
                {
                    int mod = v % (1 << bit_num);
                    v >>= bit_num;
                    temp_arr.push_back(arr[mod]);
                }
                for (int j = temp_arr.size() - 1; j >= 0 && i < _vec_len; j--, i++)
                    _vec[i].ch = temp_arr[j];
            }
        }
        else
        {
            for (int i = 0; i < _vec_len; i++)
                _vec[i].ch = arr[0];
        }
    }*/

    private static void readOtherData(BufferedReader br) {
        //read lowercase character information
       /*
        int flag;
       */

        try {
         /*           flag = Integer.parseInt(br.readLine());
            System.out.println("lowercase flag: "+flag);
            if (flag!=1){

        */
                seq_low_len = Integer.parseInt(br.readLine());
               // System.out.println("seq_low_len: "+seq_low_len);
                readPositionRangeData(br, seq_low_len , seq_low_begin, seq_low_length);
         /*   } else {
                low_loc = runLengthDecoding(br,1);//low_loc在decoding中会开辟内存，不需要在initial中开辟内存
                if(low_loc!=null)
                    seq_low_len = low_loc.length;
                else
                    seq_low_len = 0;
                diff_low_len = Integer.parseInt(br.readLine());
                System.out.println("diff_low_len: "+diff_low_len);
                readPositionRangeData(br, diff_low_len, diff_low_begin,diff_low_length);
                seqLowercaseReading(seq_low_len);
            }

          */
            //read n character information
            nCha_len = Integer.parseInt(br.readLine());
          //  System.out.println("nCha_len: "+nCha_len);
            readPositionRangeData(br, nCha_len, Main.nCha_begin, Main.nCha_length);

            //read special character information
            spe_cha_len = Integer.parseInt(br.readLine());
         //   System.out.println("spe_cha_len: "+spe_cha_len);
//            fscanf(fp, "%d", &_spe_cha_len);
            if (spe_cha_len > 0)
                readPositionRangeData(br, spe_cha_len, Main.spe_cha_pos, Main.spe_cha_ch);
//                readSpeChaData(fp, _spe_cha_len, spe_cha);
        } catch (IOException e) {
            e.printStackTrace();
        }
//        fscanf(fp, "%d", &flag);
    }

/*    private static int spaceNumber(char []str) {
        int num = 0;
        for (int i = 0;i<str.length ; i++)
            if (str[i] == ' ')
                num++;
        return num;
    }*/
    private static void readFirstMatchResult(BufferedReader br, List <MatchEntry> _mr) {
        String str;
        String []temp_str ;
        StringBuilder _misStr = new StringBuilder();
        MatchEntry _me  = new MatchEntry();
//        _mr.clear();
        int _seq_id, _pos, _length, pre_seq_id = 0, pre_pos = 0;
        try {
            while ((str=br.readLine())!=null) {
                temp_str = str.split("\\s+");
                if (temp_str.length==3) {
                    _seq_id = Integer.parseInt(temp_str[0]);
                    _pos = Integer.parseInt(temp_str[1]);
                    _length = Integer.parseInt(temp_str[2]);
//                    sscanf(temp_str, "%d%d%d", &_seq_id, &_pos, &_length);
                    _seq_id += pre_seq_id;
                    pre_seq_id = _seq_id;
                    _pos += pre_pos;
                    _length += 2;
                    pre_pos = _pos + _length;
                    getMatchResult(_seq_id, _pos, _length, _mr);
//                scanf(temp_str, "%d%d", &_pos, &_length);
                } else if (temp_str.length==2) {
                    _pos = Integer.parseInt(temp_str[0]);
                    _length = Integer.parseInt(temp_str[1]);
                    _me.setPos(_pos);
                    _me.setLength(_length);
                    _me.setMisStr(_misStr.toString());
                    _mr.add(_me);
                    _me = new MatchEntry();
                    _misStr.delete(0,_misStr.length());
                } else if (temp_str.length==1) {
                    _misStr.append(temp_str[0]);
                }
//            fgets(temp_str, 10240, fp);
            }
        }
        catch (IOException e) {
            e.printStackTrace();
        }
        //System.out.println("Read first match result complete.");
    }

    private static void getMatchResult(int _seq_id, int _pos, int _length, List <MatchEntry> _mr) {
        for (int i = 0; i < _length; i++)
            _mr.add(matchResult_vec.get(_seq_id).get(_pos++));
    }

/*    private static void readSecondMatchResult(BufferedReader br, List <MatchEntry> _mr) {
        int pre_seq_id = 0, _seq_id, pre_pos = 0, _pos, _length, spaceNum;
        char []temp_str = new char[10240];
        String _misStr = "";
        MatchEntry _me;
        _mr.clear();

        fgets(temp_str, 10240, fp);//把otherdata和matchentity之间的第一个回车抛掉
        fgets(temp_str, 10240, fp);
        try {
            _me = new MatchEntry();
            temp_str = br.readLine().toCharArray();
            //while (fgets(str, MAX_CHAR_NUM, fp) != NULL)
            while (temp_str[0] != '\n') {
                spaceNum = spaceNumber(temp_str);

                temp_str = br.readLine().toCharArray();
//            fgets(temp_str, 10240, fp);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

        System.out.println("Read second match result complete.");
    }*/

    private static void readTargetSequenceCode(List <MatchEntry> _mr) {
        int _pos, pre_pos=0, cur_pos,  _length, _seq_code_len=0, str_len;
	    char []_misStr;
        for (int i = 0; i < _mr.size();i++) {
            _pos = _mr.get(i).getPos();
            cur_pos = _pos + pre_pos;
            _length = _mr.get(i).getLength() + min_rep_len;
            pre_pos = cur_pos + _length;
            _misStr = _mr.get(i).getMisStr().toCharArray();
            //因为压缩的时候是先记录misStr，再记录matchentity，所以在解压的时候也要先恢复misStr，再恢复matchentity
            str_len = _misStr.length;
            for (int k = 0; k < str_len; k++)
                seq_code[_seq_code_len++] = integerEncoding[_misStr[k] - '0'];
            for (int m = cur_pos, n = 0; n < _length; n++, m++)
                seq_code[_seq_code_len++] = ref_code[m];
        }
//        seq_code[_seq_code_len] = '\0';
        seq_code_len = _seq_code_len;
//        System.out.println("readTargetSequenceCode complete.");
    }

    private static void saveSequenceFile(BufferedWriter wr, int _seqnum) {
        for (int i = 1; i < spe_cha_len; i++)
            spe_cha_pos[i] += spe_cha_pos[i-1];

        char []temp_seq = new char[MAX_CHA_NUM];

        int tt = 0, j = 0;
        //spe char
        for (int i = 0; i < spe_cha_len; i++) {
            //cout << spe_cha[i].pos << endl;
            while (tt < spe_cha_pos[i] && tt < seq_code_len) {
                temp_seq[j++] = seq_code[tt++];
            }
            temp_seq[j++] = (char)(spe_cha_ch[i] + 'A');
//            seq_code[j++] = (char)(spe_cha_ch[i] + 'A');
        }
        while (tt < seq_code_len) {
            temp_seq[j++] = seq_code[tt++];
        }
/*        seq_code[j] = '\0';*/
        seq_code_len = j;

        int str_len = 0;
        int r = 0;

        char []str = new char[MAX_CHA_NUM];

        //N
        for (int i = 0; i < nCha_len; i++) {
            for (j = 0; j < nCha_begin[i]; j++)
                str[str_len++] = temp_seq[r++];
            //int _length = nCha[i].length;
            for (j = 0; j < nCha_length[i]; j++)
                str[str_len++] = 'N';
        }
        while (r < seq_code_len)
            str[str_len++] = temp_seq[r++];

        //输出id
        try {
            wr.write(identifier+"\n");
        } catch (IOException e) {
            e.printStackTrace();
        }
//        fprintf(fp, "%s", _identifier);

        //lowercases
        int k = 0;
        for (int i = 0; i < seq_low_len; i++) {
            k += seq_low_begin[i];
            int temp = seq_low_length[i];
            for (j = 0; j < temp; j++) {
                str[k] = Character.toLowerCase(str[k]);
                k++;
            }
        }

        int _lineWidth = lineWidth;

        try {
//            wr.write(identifier+"\n");
            wr.flush();
//            int temp_seq_len = 0;
            for (int i = 0; i < str_len; i += _lineWidth) {
                for (j = i; j < i + _lineWidth; j++) {
                    if (j < str_len) {
                        wr.write(str[j]);
                    }
                }
                wr.write("\n");
                wr.flush();
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
//        fprintf(fp, "%s", temp_seq);
//        System.out.println("save sequence file complete.");
    }

    private static void decompress(String refPath, String namePath,String filePath ,String outPath) {
        BufferedReader br1;
        BufferedReader br2;
        BufferedWriter bw;

        readName(namePath);
        readFile(filePath);
//        readFile("J:\\gene\\output\\spark\\out2\\");D
        sec_seq_num = 45/*(int) ceil(PERCENT * seqNumber / 100)*/;
        initial();

        referenceSequenceExtraction(refPath);
        try {
            //先读行宽与id
//            br1 = new BufferedReader(new FileReader(new File(otherPath+seqName.get(seqNumber-1))));
//            readIdentifierData(br1, identifier_vec);//save identifier data, identifier_vec也是一个数组，它在initial中初始化，长度是seqNumber-1
//            lineWidth_vec = runLengthDecoding(br1,0);//read lineWidth data，decoding函数中会给lineWidth开辟内存，所以解压缩初始化时不需要给lineWidth开辟内存
            for (int i = 0; i < seqNumber; i++) {
                //读取其他信息
//                br1 = new BufferedReader(new FileReader(new File(otherPath+"\\"+seqName.get(i))));
//                br1.close();
                //读取序列信息
                br1 = new BufferedReader(new FileReader(new File(filePath+"/"+seqPath.get(i))));
                identifier = br1.readLine();
                lineWidth = Integer.parseInt(br1.readLine());
                readOtherData(br1);
                /*if (i == 0)*/ readFirstMatchResult(br1, matchResult);
//                else        readSecondMatchResult(br2, matchResult);
                if (i <= sec_seq_num && i != seqNumber - 1) matchResult_vec.add(matchResult);
                readTargetSequenceCode(matchResult);
                br1.close();

                bw = new BufferedWriter(new FileWriter(outPath+"/"+seqName.get(i)));

                saveSequenceFile(bw, i);
                matchResult = new ArrayList<>();
                bw.close();
                System.out.println("Decompressing ... （"+(i+1)+"/"+seqNumber+").");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    public static void main(String[] args) throws Exception {

        String input_ref = args[0];//参考序列路径
        String input_filePath = args[1];//已压缩文件路径(.bsc)
        String input_namePath = args[2];//源文件路径（hdfs获取文件名）
        String output_Path = args[3];//输出路径
//        String namePath = "D:\\geneEXP\\HRCMspark\\filePath.txt";
//        String filePath = "D:\\geneEXP\\HRCMspark\\output\\";
//        file = args[0];
        //只传递出压缩程序输出路径，程序将自动寻找partition文件
        long startTime = System.currentTimeMillis();
        tar.BSC(input_filePath,"/home/gene");
        decompress(input_ref,input_namePath,input_filePath,output_Path);
        System.out.println("Decompression completes. The decompression takes " + (System.currentTimeMillis() - startTime) / 1000 + "s.");
    }
}
