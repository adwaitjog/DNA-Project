#include<math.h>
#include<iostream>
#include<string>
#include<sstream>
#include<fstream>
#include<stdio.h>
#include<stdlib.h>
#include<map>
#include<vector>
#include<array>
#include<iomanip>

using namespace std;

class DNA_pool {
public:
    int add_strand(array<string,2> strand) {
        strands.push_back(strand);
    }

    string get_strand(long index) {
        return strands[index][0];
    }

    long num_strands() {
        return strands.size();
    }

private:
    vector<array<string,2>> strands;
};

//abstract encoder class with virtual encode/decode functions and helper functions implemented
class Encoder {
public:
    virtual int encode(const string& input, DNA_pool& pool) = 0;                    //abstract encode function. implement with either Goldman or XOR encoding in subclass
    virtual string decode(vector<string> strands) = 0;                    //abstract decode function; implement similarly

protected:
    short prev, primer_index = 0;                   //used to store previous nucleotide and which primer
    int payload_length = 120;                       //nucleotide length of payload
    char nucleotides[4] = {'A', 'C', 'G', 'T'};     //array of nucleotides to cycle through
    map<char,short> indices;                        //index of each nucleotide
    map<char,char> complement;
    map<short, string> ternary_map;                 //map of ASCII values to ternary string via Huffman code
    map<string, short> ascii_map;
    map<char, char> complements;
    map<string, char> xor_map;
    string primers[6][2] = {{"CGACAGTAACTACACGGCGA", "CTTGGTCAGACGAGTGCATG"}, {"GTAGCAATTGGCAGGTCCAT", "GAGTTACGCGGGGATACATG"}, {"GTAGCAATTGGCAGGTCCAT", "TGGTACGGGAACAGCACATG"}, {"CGACAGTAACTACACGGCGA", "CGTTAAGACGTAGCCCCATG"}, {"GTAGCAATTGGCAGGTCCAT", "CTCACCGCTCTTGTAGCATG"}, {"CGACAGTAACTACACGGCGA", "GACCGGCAATCTCTTCCTGG"}};                             //array of primers


    //init map of nucleotides to their indices
    int init_indices(void) {
        indices['A'] = 0;
        indices['C'] = 1;
        indices['G'] = 2;
        indices['T'] = 3;

        return 0;
    }

    //init map of nucleotides and their complements
    int init_complement(void) {
        complement['A'] = 'T';
        complement['T'] = 'A';
        complement['C'] = 'G';
        complement['G'] = 'C';

        return 0;
    }

    //init map of Huffman code
    int init_huff(void) {
        ternary_map[0] = "22201";
        ternary_map[85] = "22200";
        ternary_map[170] = "22122";
        ternary_map[127] = "22121";
        ternary_map[253] = "22120";
        ternary_map[52] = "22112";
        ternary_map[138] = "22111";
        ternary_map[41] = "22110";
        ternary_map[86] = "22102";
        ternary_map[42] = "22101";
        ternary_map[100] = "22100";
        ternary_map[44] = "22022";
        ternary_map[250] = "22020";
        ternary_map[132] = "22021";
        ternary_map[161] = "22012";
        ternary_map[98] = "22010";
        ternary_map[8] = "22002";
        ternary_map[34] = "22011";
        ternary_map[10] = "22001";
        ternary_map[149] = "22000";
        ternary_map[87] = "21222";
        ternary_map[21] = "21221";
        ternary_map[74] = "21220";
        ternary_map[36] = "21212";
        ternary_map[69] = "21210";
        ternary_map[177] = "21202";
        ternary_map[20] = "21211";
        ternary_map[213] = "21200";
        ternary_map[163] = "21201";
        ternary_map[229] = "21121";
        ternary_map[255] = "21122";
        ternary_map[197] = "21120";
        ternary_map[133] = "21112";
        ternary_map[252] = "21110";
        ternary_map[26] = "21111";
        ternary_map[173] = "21101";
        ternary_map[151] = "21102";
        ternary_map[82] = "21100";
        ternary_map[75] = "21022";
        ternary_map[37] = "21021";
        ternary_map[166] = "21011";
        ternary_map[191] = "21020";
        ternary_map[88] = "21012";
        ternary_map[63] = "21010";
        ternary_map[68] = "21001";
        ternary_map[150] = "21002";
        ternary_map[76] = "21000";
        ternary_map[4] = "20222";
        ternary_map[154] = "20221";
        ternary_map[234] = "20212";
        ternary_map[22] = "20220";
        ternary_map[162] = "20211";
        ternary_map[105] = "20210";
        ternary_map[102] = "20202";
        ternary_map[171] = "20201";
        ternary_map[104] = "20200";
        ternary_map[169] = "20122";
        ternary_map[196] = "20121";
        ternary_map[208] = "20120";
        ternary_map[84] = "20112";
        ternary_map[130] = "20111";
        ternary_map[146] = "20102";
        ternary_map[72] = "20110";
        ternary_map[16] = "20101";
        ternary_map[66] = "20100";
        ternary_map[24] = "20022";
        ternary_map[106] = "20012";
        ternary_map[223] = "20020";
        ternary_map[58] = "20021";
        ternary_map[137] = "20011";
        ternary_map[73] = "20010";
        ternary_map[101] = "20001";
        ternary_map[168] = "20002";
        ternary_map[181] = "12221";
        ternary_map[175] = "12222";
        ternary_map[251] = "20000";
        ternary_map[40] = "12220";
        ternary_map[140] = "12212";
        ternary_map[17] = "12211";
        ternary_map[83] = "12210";
        ternary_map[254] = "12202";
        ternary_map[240] = "12201";
        ternary_map[214] = "12200";
        ternary_map[53] = "12122";
        ternary_map[202] = "12112";
        ternary_map[25] = "12121";
        ternary_map[18] = "12120";
        ternary_map[247] = "12111";
        ternary_map[174] = "12110";
        ternary_map[112] = "12102";
        ternary_map[89] = "12101";
        ternary_map[210] = "12100";
        ternary_map[217] = "12012";
        ternary_map[248] = "12020";
        ternary_map[194] = "12021";
        ternary_map[182] = "12022";
        ternary_map[80] = "12011";
        ternary_map[79] = "12002";
        ternary_map[195] = "12010";
        ternary_map[12] = "12001";
        ternary_map[209] = "12000";
        ternary_map[165] = "11222";
        ternary_map[245] = "11221";
        ternary_map[2] = "11220";
        ternary_map[81] = "11212";
        ternary_map[38] = "11211";
        ternary_map[141] = "11202";
        ternary_map[211] = "11210";
        ternary_map[239] = "11200";
        ternary_map[95] = "11201";
        ternary_map[43] = "11122";
        ternary_map[224] = "11121";
        ternary_map[203] = "11112";
        ternary_map[145] = "11120";
        ternary_map[147] = "11110";
        ternary_map[19] = "11111";
        ternary_map[50] = "11101";
        ternary_map[136] = "11102";
        ternary_map[107] = "11100";
        ternary_map[134] = "11022";
        ternary_map[109] = "11021";
        ternary_map[153] = "11020";
        ternary_map[148] = "11002";
        ternary_map[205] = "11010";
        ternary_map[212] = "11011";
        ternary_map[54] = "11012";
        ternary_map[241] = "11000";
        ternary_map[156] = "11001";
        ternary_map[115] = "10222";
        ternary_map[116] = "10221";
        ternary_map[78] = "10220";
        ternary_map[67] = "10211";
        ternary_map[70] = "10212";
        ternary_map[178] = "10210";
        ternary_map[159] = "10202";
        ternary_map[142] = "10201";
        ternary_map[92] = "10200";
        ternary_map[48] = "10122";
        ternary_map[90] = "10120";
        ternary_map[218] = "10121";
        ternary_map[126] = "10112";
        ternary_map[39] = "10111";
        ternary_map[219] = "10102";
        ternary_map[167] = "10110";
        ternary_map[114] = "10101";
        ternary_map[172] = "10022";
        ternary_map[14] = "10100";
        ternary_map[120] = "10020";
        ternary_map[139] = "10021";
        ternary_map[160] = "10012";
        ternary_map[33] = "10011";
        ternary_map[179] = "10010";
        ternary_map[117] = "10002";
        ternary_map[225] = "10001";
        ternary_map[129] = "10000";
        ternary_map[183] = "02222";
        ternary_map[230] = "02220";
        ternary_map[35] = "02221";
        ternary_map[93] = "02210";
        ternary_map[6] = "02211";
        ternary_map[32] = "02212";
        ternary_map[56] = "02201";
        ternary_map[158] = "02202";
        ternary_map[185] = "02121";
        ternary_map[47] = "02122";
        ternary_map[143] = "02200";
        ternary_map[123] = "02111";
        ternary_map[204] = "02120";
        ternary_map[242] = "02112";
        ternary_map[111] = "02110";
        ternary_map[103] = "02102";
        ternary_map[108] = "02101";
        ternary_map[9] = "02100";
        ternary_map[65] = "02022";
        ternary_map[249] = "02020";
        ternary_map[13] = "02021";
        ternary_map[180] = "02012";
        ternary_map[226] = "02001";
        ternary_map[144] = "02002";
        ternary_map[15] = "02010";
        ternary_map[57] = "02011";
        ternary_map[128] = "02000";
        ternary_map[135] = "01220";
        ternary_map[243] = "01221";
        ternary_map[190] = "01222";
        ternary_map[207] = "01212";
        ternary_map[77] = "01211";
        ternary_map[45] = "01210";
        ternary_map[91] = "01202";
        ternary_map[192] = "01201";
        ternary_map[186] = "01122";
        ternary_map[216] = "01200";
        ternary_map[97] = "01112";
        ternary_map[118] = "01120";
        ternary_map[246] = "01121";
        ternary_map[215] = "01111";
        ternary_map[51] = "01102";
        ternary_map[206] = "01110";
        ternary_map[184] = "01100";
        ternary_map[227] = "01101";
        ternary_map[233] = "01022";
        ternary_map[237] = "01021";
        ternary_map[188] = "01020";
        ternary_map[113] = "01012";
        ternary_map[49] = "01011";
        ternary_map[201] = "01010";
        ternary_map[155] = "01002";
        ternary_map[222] = "01000";
        ternary_map[231] = "01001";
        ternary_map[5] = "00222";
        ternary_map[27] = "00221";
        ternary_map[131] = "00212";
        ternary_map[164] = "00220";
        ternary_map[3] = "00211";
        ternary_map[46] = "00210";
        ternary_map[119] = "00201";
        ternary_map[28] = "00202";
        ternary_map[176] = "00200";
        ternary_map[23] = "00122";
        ternary_map[64] = "00121";
        ternary_map[157] = "00120";
        ternary_map[187] = "00112";
        ternary_map[244] = "00110";
        ternary_map[238] = "00111";
        ternary_map[96] = "00102";
        ternary_map[235] = "00101";
        ternary_map[60] = "00022";
        ternary_map[1] = "00100";
        ternary_map[110] = "00021";
        ternary_map[200] = "00011";
        ternary_map[221] = "00020";
        ternary_map[99] = "00012";
        ternary_map[31] = "00010";
        ternary_map[198] = "00002";
        ternary_map[193] = "00001";
        ternary_map[125] = "00000";
        ternary_map[124] = "222222";
        ternary_map[152] = "222221";
        ternary_map[122] = "222220";
        ternary_map[71] = "222212";
        ternary_map[94] = "222211";
        ternary_map[220] = "222210";
        ternary_map[29] = "222202";
        ternary_map[199] = "222201";
        ternary_map[61] = "222200";
        ternary_map[11] = "222122";
        ternary_map[228] = "222121";
        ternary_map[62] = "222120";
        ternary_map[55] = "222112";
        ternary_map[121] = "222111";
        ternary_map[7] = "222110";
        ternary_map[30] = "222102";
        ternary_map[232] = "222101";
        ternary_map[189] = "222100";
        ternary_map[59] = "222021";
        ternary_map[236] = "222022";

        return 0;
    }

    int init_ascii(void) {
        ascii_map["22201"] = 0;
        ascii_map["22200"] = 85;
        ascii_map["22122"] = 170;
        ascii_map["22121"] = 127;
        ascii_map["22120"] = 253;
        ascii_map["22112"] = 52;
        ascii_map["22111"] = 138;
        ascii_map["22110"] = 41;
        ascii_map["22102"] = 86;
        ascii_map["22101"] = 42;
        ascii_map["22100"] = 100;
        ascii_map["22022"] = 44;
        ascii_map["22020"] = 250;
        ascii_map["22021"] = 132;
        ascii_map["22012"] = 161;
        ascii_map["22010"] = 98;
        ascii_map["22002"] = 8;
        ascii_map["22011"] = 34;
        ascii_map["22001"] = 10;
        ascii_map["22000"] = 149;
        ascii_map["21222"] = 87;
        ascii_map["21221"] = 21;
        ascii_map["21220"] = 74;
        ascii_map["21212"] = 36;
        ascii_map["21210"] = 69;
        ascii_map["21202"] = 177;
        ascii_map["21211"] = 20;
        ascii_map["21200"] = 213;
        ascii_map["21201"] = 163;
        ascii_map["21121"] = 229;
        ascii_map["21122"] = 255;
        ascii_map["21120"] = 197;
        ascii_map["21112"] = 133;
        ascii_map["21110"] = 252;
        ascii_map["21111"] = 26;
        ascii_map["21101"] = 173;
        ascii_map["21102"] = 151;
        ascii_map["21100"] = 82;
        ascii_map["21022"] = 75;
        ascii_map["21021"] = 37;
        ascii_map["21011"] = 166;
        ascii_map["21020"] = 191;
        ascii_map["21012"] = 88;
        ascii_map["21010"] = 63;
        ascii_map["21001"] = 68;
        ascii_map["21002"] = 150;
        ascii_map["21000"] = 76;
        ascii_map["20222"] = 4;
        ascii_map["20221"] = 154;
        ascii_map["20212"] = 234;
        ascii_map["20220"] = 22;
        ascii_map["20211"] = 162;
        ascii_map["20210"] = 105;
        ascii_map["20202"] = 102;
        ascii_map["20201"] = 171;
        ascii_map["20200"] = 104;
        ascii_map["20122"] = 169;
        ascii_map["20121"] = 196;
        ascii_map["20120"] = 208;
        ascii_map["20112"] = 84;
        ascii_map["20111"] = 130;
        ascii_map["20102"] = 146;
        ascii_map["20110"] = 72;
        ascii_map["20101"] = 16;
        ascii_map["20100"] = 66;
        ascii_map["20022"] = 24;
        ascii_map["20012"] = 106;
        ascii_map["20020"] = 223;
        ascii_map["20021"] = 58;
        ascii_map["20011"] = 137;
        ascii_map["20010"] = 73;
        ascii_map["20001"] = 101;
        ascii_map["20002"] = 168;
        ascii_map["12221"] = 181;
        ascii_map["12222"] = 175;
        ascii_map["20000"] = 251;
        ascii_map["12220"] = 40;
        ascii_map["12212"] = 140;
        ascii_map["12211"] = 17;
        ascii_map["12210"] = 83;
        ascii_map["12202"] = 254;
        ascii_map["12201"] = 240;
        ascii_map["12200"] = 214;
        ascii_map["12122"] = 53;
        ascii_map["12112"] = 202;
        ascii_map["12121"] = 25;
        ascii_map["12120"] = 18;
        ascii_map["12111"] = 247;
        ascii_map["12110"] = 174;
        ascii_map["12102"] = 112;
        ascii_map["12101"] = 89;
        ascii_map["12100"] = 210;
        ascii_map["12012"] = 217;
        ascii_map["12020"] = 248;
        ascii_map["12021"] = 194;
        ascii_map["12022"] = 182;
        ascii_map["12011"] = 80;
        ascii_map["12002"] = 79;
        ascii_map["12010"] = 195;
        ascii_map["12001"] = 12;
        ascii_map["12000"] = 209;
        ascii_map["11222"] = 165;
        ascii_map["11221"] = 245;
        ascii_map["11220"] = 2;
        ascii_map["11212"] = 81;
        ascii_map["11211"] = 38;
        ascii_map["11202"] = 141;
        ascii_map["11210"] = 211;
        ascii_map["11200"] = 239;
        ascii_map["11201"] = 95;
        ascii_map["11122"] = 43;
        ascii_map["11121"] = 224;
        ascii_map["11112"] = 203;
        ascii_map["11120"] = 145;
        ascii_map["11110"] = 147;
        ascii_map["11111"] = 19;
        ascii_map["11101"] = 50;
        ascii_map["11102"] = 136;
        ascii_map["11100"] = 107;
        ascii_map["11022"] = 134;
        ascii_map["11021"] = 109;
        ascii_map["11020"] = 153;
        ascii_map["11002"] = 148;
        ascii_map["11010"] = 205;
        ascii_map["11011"] = 212;
        ascii_map["11012"] = 54;
        ascii_map["11000"] = 241;
        ascii_map["11001"] = 156;
        ascii_map["10222"] = 115;
        ascii_map["10221"] = 116;
        ascii_map["10220"] = 78;
        ascii_map["10211"] = 67;
        ascii_map["10212"] = 70;
        ascii_map["10210"] = 178;
        ascii_map["10202"] = 159;
        ascii_map["10201"] = 142;
        ascii_map["10200"] = 92;
        ascii_map["10122"] = 48;
        ascii_map["10120"] = 90;
        ascii_map["10121"] = 218;
        ascii_map["10112"] = 126;
        ascii_map["10111"] = 39;
        ascii_map["10102"] = 219;
        ascii_map["10110"] = 167;
        ascii_map["10101"] = 114;
        ascii_map["10022"] = 172;
        ascii_map["10100"] = 14;
        ascii_map["10020"] = 120;
        ascii_map["10021"] = 139;
        ascii_map["10012"] = 160;
        ascii_map["10011"] = 33;
        ascii_map["10010"] = 179;
        ascii_map["10002"] = 117;
        ascii_map["10001"] = 225;
        ascii_map["10000"] = 129;
        ascii_map["02222"] = 183;
        ascii_map["02220"] = 230;
        ascii_map["02221"] = 35;
        ascii_map["02210"] = 93;
        ascii_map["02211"] = 6;
        ascii_map["02212"] = 32;
        ascii_map["02201"] = 56;
        ascii_map["02202"] = 158;
        ascii_map["02121"] = 185;
        ascii_map["02122"] = 47;
        ascii_map["02200"] = 143;
        ascii_map["02111"] = 123;
        ascii_map["02120"] = 204;
        ascii_map["02112"] = 242;
        ascii_map["02110"] = 111;
        ascii_map["02102"] = 103;
        ascii_map["02101"] = 108;
        ascii_map["02100"] = 9;
        ascii_map["02022"] = 65;
        ascii_map["02020"] = 249;
        ascii_map["02021"] = 13;
        ascii_map["02012"] = 180;
        ascii_map["02001"] = 226;
        ascii_map["02002"] = 144;
        ascii_map["02010"] = 15;
        ascii_map["02011"] = 57;
        ascii_map["02000"] = 128;
        ascii_map["01220"] = 135;
        ascii_map["01221"] = 243;
        ascii_map["01222"] = 190;
        ascii_map["01212"] = 207;
        ascii_map["01211"] = 77;
        ascii_map["01210"] = 45;
        ascii_map["01202"] = 91;
        ascii_map["01201"] = 192;
        ascii_map["01122"] = 186;
        ascii_map["01200"] = 216;
        ascii_map["01112"] = 97;
        ascii_map["01120"] = 118;
        ascii_map["01121"] = 246;
        ascii_map["01111"] = 215;
        ascii_map["01102"] = 51;
        ascii_map["01110"] = 206;
        ascii_map["01100"] = 184;
        ascii_map["01101"] = 227;
        ascii_map["01022"] = 233;
        ascii_map["01021"] = 237;
        ascii_map["01020"] = 188;
        ascii_map["01012"] = 113;
        ascii_map["01011"] = 49;
        ascii_map["01010"] = 201;
        ascii_map["01002"] = 155;
        ascii_map["01000"] = 222;
        ascii_map["01001"] = 231;
        ascii_map["00222"] = 5;
        ascii_map["00221"] = 27;
        ascii_map["00212"] = 131;
        ascii_map["00220"] = 164;
        ascii_map["00211"] = 3;
        ascii_map["00210"] = 46;
        ascii_map["00201"] = 119;
        ascii_map["00202"] = 28;
        ascii_map["00200"] = 176;
        ascii_map["00122"] = 23;
        ascii_map["00121"] = 64;
        ascii_map["00120"] = 157;
        ascii_map["00112"] = 187;
        ascii_map["00110"] = 244;
        ascii_map["00111"] = 238;
        ascii_map["00102"] = 96;
        ascii_map["00101"] = 235;
        ascii_map["00022"] = 60;
        ascii_map["00100"] = 1;
        ascii_map["00021"] = 110;
        ascii_map["00011"] = 200;
        ascii_map["00020"] = 221;
        ascii_map["00012"] = 99;
        ascii_map["00010"] = 31;
        ascii_map["00002"] = 198;
        ascii_map["00001"] = 193;
        ascii_map["00000"] = 125;
        ascii_map["222222"] = 124;
        ascii_map["222221"] = 152;
        ascii_map["222220"] = 122;
        ascii_map["222212"] = 71;
        ascii_map["222211"] = 94;
        ascii_map["222210"] = 220;
        ascii_map["222202"] = 29;
        ascii_map["222201"] = 199;
        ascii_map["222200"] = 61;
        ascii_map["222122"] = 11;
        ascii_map["222121"] = 228;
        ascii_map["222120"] = 62;
        ascii_map["222112"] = 55;
        ascii_map["222111"] = 121;
        ascii_map["222110"] = 7;
        ascii_map["222102"] = 30;
        ascii_map["222101"] = 232;
        ascii_map["222100"] = 189;
        ascii_map["222021"] = 59;
        ascii_map["222022"] = 236;

        return 0;
    }

    int init_primers(void) {
        array<array<string,2>, 6> primers = {{
            {"CGACAGTAACTACACGGCGA", "CTTGGTCAGACGAGTGCATG"},
            {"GTAGCAATTGGCAGGTCCAT", "GAGTTACGCGGGGATACATG"},
            {"GTAGCAATTGGCAGGTCCAT", "TGGTACGGGAACAGCACATG"},
            {"CGACAGTAACTACACGGCGA", "CGTTAAGACGTAGCCCCATG"},
            {"GTAGCAATTGGCAGGTCCAT", "CTCACCGCTCTTGTAGCATG"},
            {"CGACAGTAACTACACGGCGA", "GACCGGCAATCTCTTCCTGG"},
        }};
        return 0;
    }

    int init_xor_map() {
        xor_map["AxA"] = 'A';
        xor_map["AxC"] = 'C';
        xor_map["CxA"] = 'C';
        xor_map["AxG"] = 'G';
        xor_map["GxA"] = 'G';
        xor_map["AxT"] = 'T';
        xor_map["TxA"] = 'T';
        xor_map["CxC"] = 'A';
        xor_map["CxG"] = 'T';
        xor_map["GxC"] = 'T';
        xor_map["CxT"] = 'G';
        xor_map["TxC"] = 'G';
        xor_map["GxG"] = 'A';
        xor_map["GxT"] = 'C';
        xor_map["TxG"] = 'C';
        xor_map["TxT"] = 'A';

        return 0;
    }

    int init_stuff(void) {
        init_indices();
        init_huff();
        init_ascii();
        init_complement();
        init_primers();
        init_xor_map();
        return 1;
    }

    string to_string(int num) {
        ostringstream convert;
        convert << num;
        return convert.str();
    }

    template <class BidirectionalIterator> void reverse (BidirectionalIterator first, BidirectionalIterator last)
    {
        while ((first!=last)&&(first!=--last)) {
            std::iter_swap (first,last);
            ++first;
        }
    }

    /* Simple base 10 to base X converter with input parameters (X, base).
    This prints at at a digit limit of 4, but can be modified. This means your
    input can be no larger than base^{digit_limit} */
    string base_converter (int input, int base) {
        string output;
        div_t divresult;
        int remainder; int digit_limit = 4;
        for (int i = 0; i < digit_limit; i++) {
            if (i == 0) { divresult = div(input, base); }
            else {
                divresult = div(divresult.quot, base);
            }
            remainder = divresult.rem; //sets remainder in ternary form
            output.append(to_string(remainder));
        }
        reverse(output.begin(), output.end()); //must reverse order of remainders
        return output;
    }

    //ternary
    string to_ternary(const string& input) {
        string ternary = "";

        //create string of ternary digits
        for(int i = 0; i < input.length(); i++)
            ternary += ternary_map[input[i]];

        return ternary;
    }

    //rotating encoding of nucleotides
    char rotate_nucleotides(short t) {
        prev = (prev+1+t)%4;
        return nucleotides[prev];
    }

    //read in string, convert each ternary digit to a nucleotide
    //return nucleotide sequence representing input stream
    string to_nucleotides(const string& ternary) {
        string n = "";

        //create nucleotide payload
        for(int j = 0; j < ternary.length(); j++)
            n += rotate_nucleotides(ternary[j]-48);

        return n;
    }

    string get_address(int num) {
        string ternary = base_converter(num, 3);
        string address = string(6-ternary.length(), 0) + ternary;
        return to_nucleotides(ternary);
    }

    string to_complement(const string& strand) {
        string antisense;
        for(int i = 0; i < strand.length(); i++) {
            antisense += complement[strand[i]];
        }
        return antisense;
    }

    //generate one strand of DNA
    map<string,string> synthesize_strand(int is_complement, const string& ternary, int num) {
        map<string,string> strand;
        prev = indices[primers[primer_index][0].back()];
        strand["s"] = rotate_nucleotides(is_complement);
        strand["payload"] = to_nucleotides(ternary);
        strand["address"] = get_address(num);

        return strand;
    }

    //store information in DNA
    vector<map<string,string>> to_dna(const string& input) {
        string ternary;

        vector<map<string,string>> ret_list;

        ternary = to_ternary(input);

        short num_strands = floor(ternary.length()/payload_length)+1;
        short strand_size = ceil(ternary.length()/num_strands);

        for (int i = 0; i < num_strands-1; i++)
            ret_list.push_back(synthesize_strand(0, ternary.substr(i*strand_size, strand_size), i));
        ret_list.push_back(synthesize_strand(0, ternary.substr((num_strands-1)*strand_size), num_strands-1));

        return ret_list;
    }

    int to_pool(vector<map<string,string>> strands, DNA_pool& pool) {
        for(map<string,string> s : strands) {
            array<string,2> strand;
            strand[0] = primers[primer_index][0] + s["s"] + s["payload"] + s["address"] + s["s"] + primers[primer_index][1];
            strand[1] = to_complement(strand[0]);
            pool.add_strand(strand);
        }
        primer_index++;
        return 1;
    }

    const string from_ternary(const string& ternary_string) {
        string decoded_payload;
        int i = 0;
        string byte;
        //loop till we've found all bytes
        while (i < ternary_string.length()) {
            //check if byte is 5 digits long, if not, find its 6 digit makeup
            byte = ternary_string.substr(i, 5);
            if (ascii_map.find(byte) != ascii_map.end() ) { i+=5; } // m.find("f") == m.end()
            else {
                byte = ternary_string.substr(i, 6);
                i+=6;
            }
            //cout << byte + "\t";
            decoded_payload += (char)ascii_map[byte];
            if (i > ternary_string.length()) { cout << "ERROR   \n";}
        }
        return decoded_payload;
    }

    //subfunction, finds index of a value in an array
    int get_index(char *ra, size_t len, char c) {
      char* last = ra + len;
      //char* get = find(ra, last, c);
      char* get = ra;
      return (last == get)? -1 : (get - ra);
    }

    //reads the last nucleotide and the nucleotide in question to determine
    //its appropriate ternary bit equivalent
    int sequence_nucleotide(char c) {
        int index = indices[c];
        int ternary = (index - prev + 3) % 4;
        prev = index;
        //std::cout << ternary << std::endl;
        return ternary;
    }

    //calls sequence nucleotide, forloops each N and returns ternary bits
    string from_nucleotides(const string& N) {
        string ternary_string;
        //cout << endl << "N     " << N << endl;
        //iterates through payload, converting nucleotide array
        //into ternary array
        for (int i = 0; i < N.length(); i++) {
            ostringstream convert;
            convert << sequence_nucleotide(N[i]);
            ternary_string += convert.str();
        }
        return ternary_string;
    }

    string sequence(const string& input) {
        string output;
        //prev = 0; //reset prev
        //cout << endl << "string to be decoded: " << input << endl;
        output = from_nucleotides(input);
        //cout << output;
        //cout << endl << "from_nucleotides: " << output << endl;
        return output;
    }
};

//Encoder for testing purposes once we're ready to test object-oriented approach
class TestingEncoder : Encoder {
public:
    TestingEncoder () {
        init_stuff();
    }

    int encode(const string& input, DNA_pool& pool) {
        to_pool(to_dna(input), pool);
        return 1;
    }

    string decode(vector<string> strands) {
        string strand, ternary = "";
        for (int i = 0; i < strands.size(); i++) {
            strand = strands[i];
            prev = indices[strand[20]];
            ternary += sequence(strand.substr(21, strand.length() - 46));
        }
        return from_ternary(ternary);
    }
};

/* Goldman encoder/decoder class
class GoldmanEncoder : Encoder {
public:
    //add
}; */

 //XOR encoder/decoder class
class XOREncoder : Encoder {
public:
    XOREncoder () {
        init_stuff();
    }

    string dna_xor(string payload1, string payload2) {
        string new_payload = "";
        for(int i = 0; i < payload1.length() && i < payload2.length(); i++) {
            new_payload += xor_map[payload1.substr(i,1) + "x" + payload2.substr(i,1)];
        }
        return new_payload;
    }

    string xor_address(int num1, int num2) {
        //return to_nucleotides("1" + base_converter(num1,3).substr(5) + base_converter(num2,3).substr(5));
        return to_nucleotides("1" + get_address(num1).substr(1));
    }

    vector<map<string,string>> xor_encode(vector<map<string,string>> strands) {
        int len = strands.size();
        for(int i = 0; i+1 < len; i+=2) {
            map<string,string> xor_strand;
            prev = indices[primers[primer_index][0].back()];
            xor_strand["s"] = rotate_nucleotides(0);
            xor_strand["payload"] = dna_xor(strands[i]["payload"], strands[i+1]["payload"]);
            xor_strand["address"] = xor_address(i, i+1);
            strands.push_back(xor_strand);
        }
        return strands;
    }

    int encode(const string& input, DNA_pool& pool) {
        to_pool(xor_encode(to_dna(input)), pool);
        return 1;
    }

    string decode(vector<string> strands) {
        string strand, ternary = "";
        int num_strands = ceil(.6*strands.size());
        for (int i = 0; i < num_strands; i++) {
            strand = strands[i];
            prev = indices[strand[20]];
            ternary += sequence(strand.substr(21, strand.length() - 46));
        }
        return from_ternary(ternary);
    }
};

int main(void) {
    //string test = "single strand test";
    string test = "This is a paragraph of text that I'm using to test DNA synthesis simulation. This should result in multiple strands being synthesized to test that functionality. Testing 123...";
    //string test = "Happy birthday Halle!";
    //TestingEncoder *tester = new TestingEncoder;
    XOREncoder *tester = new XOREncoder;
    DNA_pool *pool = new DNA_pool;
    tester->encode(test, *pool);
    string out_string, strand, ret_string;
    vector<string> strands;
    cout << "In nucleotides: \n";
    for(int i = 0; i < pool->num_strands(); i++) {
        strand = pool->get_strand(i);
        out_string = strand + "\n\n";
        cout << out_string;
        strands.push_back(strand);
    }
    cout << "Original data: \n";
    ret_string = tester->decode(strands);
    cout << ret_string;
    return 0;
}
