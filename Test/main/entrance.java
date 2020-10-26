package main;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

public class entrance {

    private static int MAX_SEQ_NUM = 2000;//maximum sequence number
    private static int MAX_CHA_NUM = 1 << 28;//maximum length of a chromosome
    private static int LINE_CHA_NUM =200;
    private static int PERCENT = 10; //the percentage of compressed sequence number uses as reference
    private static int kMerLen = 12; //the length of k-mer
    private static int kmer_bit_num = 2 * kMerLen; //bit numbers of k-mer
    private static int hashTableLen = 1 << kmer_bit_num; // length of hash table
    private static int VEC_SIZE = 1 <<20; //length for other character arrays
    private static int min_rep_len = 25;   //minimum replace length, matched string length exceeds min_rep_len, saved as matched information
    private static long startTime;

    private static String identifier;
    private static int lineWidth, ref_code_len, seq_code_len, ref_low_len, seq_low_len, diff_low_len, nCha_len, spe_cha_len, seqNumber, seqBucketLen;
    private static int sec_seq_num; //the referenced sequence number used for second compress
    private static char ref_code[]=new char[MAX_CHA_NUM];
    private static char seq_code[]= new char[MAX_CHA_NUM];//mismatched subsequence
/*    private static int refLoc[]= new int[MAX_CHA_NUM]; //reference hash location
    private static int refBucket[]= new int[hashTableLen]; //reference hash bucket*/
    private static int low_loc[]= new int[VEC_SIZE]; //lowercase tuple location

    private static int ref_low_begin[];
    private static int ref_low_length[];
    private static int diff_low_begin[];
    private static int diff_low_length[];
    private static int nCha_begin[];
    private static int nCha_length[];
    private static int spe_cha_pos[];
    private static int spe_cha_ch[];
    private static int seq_low_begin[];
    private static int seq_low_length[];
    private static int lineWidth_vec[] = new int[2000];

    private static List<String> SeqName;
    private static List<String> identifier_vec;

    private static void initial() {
        ref_low_begin = new int[VEC_SIZE];
        ref_low_length = new int[VEC_SIZE];
        diff_low_begin = new int[VEC_SIZE];
        diff_low_length = new int[VEC_SIZE];
        nCha_begin = new int[VEC_SIZE/2];
        nCha_length = new int[VEC_SIZE/2];
        spe_cha_ch = new int[VEC_SIZE/2];
        spe_cha_pos = new int[VEC_SIZE/2];
        seq_low_begin = new int[VEC_SIZE];
        seq_low_length = new int[VEC_SIZE];

        identifier_vec = new ArrayList<>();
        SeqName = new ArrayList<>();
        //lineWidth_vec.reverse();
        //seqBucket_vec.reserve(seqNumber);
        //seqLoc_vec.reserve(seqNumber);
    }

    private static void referenceSequenceExtraction(String str_referenceName) {
        int _seq_code_len = 0, _ref_low_len = 1, letters_len = 0;//record lowercase from 1, diff_lowercase_loc[i]=0 means mismatching
        char temp_cha;
        boolean flag = true;
        String str;
        char cha[];      //the content of one line
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
        System.out.println("Extraction of reference sequence complete. Reference code length: %d. Lowercase length: %d.\n"+ ref_code_len+"\t"+ref_low_len);

    }

    private static void seqLowercaseReading(int _seq_low_len, int _diff_low_len) {
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

    private static int[] runLengthDecoding(BufferedReader br, int tolerance) //这个地方是*&vec，但是原来都是**vec，后面看下会出不会出问题
    {
//        File fp = new File(filename);
        int code_len, temp;
        int[] vec;
        List<Integer> code = new ArrayList<>();
        try {
//            BufferedReader br = new BufferedReader(new FileReader(fp));
            code_len = Integer.parseInt(br.readLine());
            for (int i = 0; i < code_len; i++) {
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
            code_len = Integer.parseInt(br.readLine());//读取identifier长度
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

    private static void readOtherData(BufferedReader br) {
        //read lowercase character information
        int flag;
        try {
            flag = Integer.parseInt(br.readLine());
            if (flag!=1){
                seq_low_len = Integer.parseInt(br.readLine());
                if(seq_low_len!=0)
                    readPositionRangeData(br, seq_low_len , seq_low_begin, seq_low_length);
            } else {
                low_loc = runLengthDecoding(br,1);//low_loc在decoding中会开辟内存，不需要在initial中开辟内存
                readPositionRangeData(br, diff_low_len, diff_low_begin, diff_low_length);
                seqLowercaseReading(seq_low_len, diff_low_len);
            }
            //read n character information
            nCha_len = Integer.parseInt(br.readLine());
            readPositionRangeData(br, nCha_len, nCha_begin, nCha_length);

            //read special character information
            spe_cha_len = Integer.parseInt(br.readLine());
//            fscanf(fp, "%d", &_spe_cha_len);
            if (spe_cha_len > 0)
                readPositionRangeData(br, spe_cha_len, spe_cha_ch, spe_cha_pos);
//                readSpeChaData(fp, _spe_cha_len, spe_cha);
        } catch (IOException e) {
            e.printStackTrace();
        }
//        fscanf(fp, "%d", &flag);
    }

    private static void readFile(String filePath) throws NullPointerException{
        File fp =new File(filePath);
        String []names = fp.list();
        for (String a: names) {
            if(!a.contains(".crc")&&a.contains("part-0"))
                SeqName.add(a);
        }
//        seqNumber = SeqName.size()-1;
/*        String temp_name;
        //the length of filename
        BufferedReader br;
        try {
            br=new BufferedReader(new InputStreamReader(new DataInputStream(new FileInputStream(fp))));
            while ((temp_name=br.readLine())!=null) {
                seqName.add(temp_name);//放在对象末尾
            }
            br.close();

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }*/
        seqNumber = SeqName.size();
    }


    public static void main(String[] args) {
        initial();
        String otherPath = "J:\\gene\\output\\spark\\out1";
        String filePath = "J:\\gene\\output\\spark\\out2";
        readFile(otherPath);
        referenceSequenceExtraction("J:\\gene\\hg13\\chr1.fa");
        BufferedReader br1;
        BufferedReader br2;
        try {
            for (int i = 0; i < seqNumber-1; i++){
                br1 = new BufferedReader(new FileReader(new File(otherPath+SeqName.get(i))));
                readOtherData(br1);
                br1.close();
                br2 = new BufferedReader(new FileReader(new File(filePath+SeqName.get(i))));

/*                if (i == 0) readFirstMatchResult(br2, matchResult);
                else        readSecondMatchResult(br2, matchResult);
                if (i <= sec_seq_num && i != seqNumber - 2) matchResult_vec.add(matchResult);
                readTargetSequenceCode(matchResult);*/

            }
            br1 = new BufferedReader(new FileReader(new File(otherPath+SeqName.get(seqNumber-1))));
            lineWidth_vec = runLengthDecoding(br1,0);//read lineWidth data，decoding函数中会给lineWidth开辟内存，所以解压缩初始化时不需要给lineWidth开辟内存
            readIdentifierData(br1, identifier_vec);//save identifier data, identifier_vec也是一个数组，它在initial中初始化，长度是seqNumber-1
            br1.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
