using System.Collections;
using System.Collections.Generic;
using System.Linq;
using UnityEngine;

public class Main : MonoBehaviour {

    //定数
    public static double PCL_DST = 0.02;

    public static double MIN_X = (0.0 - PCL_DST * 3);
    public static double MIN_Y = (0.0 - PCL_DST * 3);
    public static double MIN_Z = (0.0 - PCL_DST * 3);

    public static double MAX_X = (1.0 + PCL_DST * 3);
    public static double MAX_Y = (0.2 + PCL_DST * 3);
    public static double MAX_Z = (0.6 + PCL_DST * 30);

    public static int GST = -1;
    public static int FLD = 0;
    public static int WLL = 1;
    public static int NUM_TYP = 2;

    public static int DNS_FLD = 1000;
    public static int DNS_WLL = 1000;

    public static double DT = 0.0005;
    public static double FIN_TIM = 1.0;

    public static double SND = 22.0;
    public static double OPT_FQC = 100;
    public static double KNM_VSC = 0.000001;

    public static int DIM = 3;

    public static double CRT_NUM = 0.1;
    public static double COL_RAT = 0.2;
    public static double DST_LMT_RAT = 0.9;

    public static double G_X = 0.0;
    public static double G_Y = 0.0;
    public static double G_Z = -9.8;

    public double wei(double dist,double re)
    {
        return ((re / dist) - 1.0);
    }

    public static int iLP;
    public static int iF;
    public static double TIM;
    public static int nP;

    public static List<Vector3> Acc = new List<Vector3>();
    public static List<Vector3> Pos = new List<Vector3>();
    public static List<Vector3> Vel = new List<Vector3>();
    public static List<double> Prs = new List<double>();
    public static List<double> pav = new List<double>();
    public static List<int> Typ = new List<int>();

    public static List<int> bfst = new List<int>(nBxyz);//バケットに格納された先頭の粒子番号 バケット数分の配列
    public static List<int> blst = new List<int>(nBxyz);//バケットに格納された最後尾の粒子番号 バケット数分の配列
    public static List<int> nxt = new List<int>(nP);//同じバケット内の次の粒子番号 粒子数分の配列


    public static double r;
    public static double r2;

    public static double DB;
    public static double DB2;
    public static double DBinv;

    public static int nBx;
    public static int nBy;
    public static int nBz;
    public static int nBxy;
    public static int nBxyz;



    public static double n0;
    public static double lmd;
    public static double A1;
    public static double A2;
    public static double A3;

    public static double rlim;
    public static double rlim2;
    public static double COL;

    public static List<double> Dns = new List<double>();
    public static List<double> invDns = new List<double>();

    


    List<GameObject> list_particle_ = new List<GameObject>();

    void RdDat(){

        Debug.Log("read input");
        TextAsset textasset = new TextAsset(); //テキストファイルのデータを取得するインスタンスを作成
        textasset = Resources.Load("dambreak", typeof(TextAsset)) as TextAsset; //Resourcesフォルダから対象テキストを取得
        string TextLines = textasset.text; //テキスト全体をstring型で入れる変数を用意して入れる
        string[] textMessage; //テキストの加工前の一行を入れる変数
                                      //Splitで一行づつを代入した1次配列を作成
        textMessage = TextLines.Split('\n'); //

        float px = 0;
        float py = 0;
        float pz = 0;

        double vx = 0;
        double vy = 0;
        double vz = 0;

        double lx = 0;
        double ly = 0;
        double lz = 0;

        //行数を取得
        int rowLength = textMessage.Length;
        Debug.Log(textMessage[0]);
        nP = int.Parse(textMessage[0]);
        for (int i = 0; i < nP+1; i++)
        {
            if(i == 0){
                
            }
            else
            {
                Typ.Add(int.Parse(textMessage[i].Split(' ')[2]));
                px = (float)(double.Parse(textMessage[i].Split(' ')[3]));
                py = (float)(double.Parse(textMessage[i].Split(' ')[4]));
                pz = (float)(double.Parse(textMessage[i].Split(' ')[5]));
                Vector3 p = new Vector3(px, py, pz);
                Pos.Add(p);
                vx = (float)(double.Parse(textMessage[i].Split(' ')[6]));
                vy = (float)(double.Parse(textMessage[i].Split(' ')[7]));
                vz = (float)(double.Parse(textMessage[i].Split(' ')[8]));
                Vector3 v = new Vector3((float)vx, (float)vy, (float)vz);
                Vel.Add(v);
                Prs.Add((float)(double.Parse(textMessage[i].Split(' ')[9])));
                pav.Add((float)(double.Parse(textMessage[i].Split(' ')[10])));
                Vector3 a = new Vector3(0.0f, 0.0f, 0.0f);
                Acc.Add(a);

                
                

                
                
            }
            

        }
        Debug.Log("Pos[].count " + Pos.Count());

    }
    void SetPara()
    {
        double tn0 = 0.0;
        double tlmd = 0.0;
        //粒子を仮想的に配置
        for(int ix = -4;ix < 5; ix++)
        {
            for (int iy = -4; ix < 5; ix++)
            {
                for (int iz = -4; ix < 5; ix++)
                {
                    double x = PCL_DST * ix;
                    double y = PCL_DST * iy;
                    double z = PCL_DST * iz;
                    double dist2 = x * x + y * y + z * z;
                    if(dist2 <= r2)
                    {
                        if (dist2 == 0.0) continue;
                        double dist = Mathf.Sqrt((float)dist2);
                        tn0 += wei(dist, r);
                        tlmd += dist2 * wei(dist, r);                        
                    }

                }
            }
        }
        n0 = tn0;           //初期粒子数密度
        lmd = tlmd / tn0;   //ラプラシアンモデルの係数λ
        A1 = 2.0 * KNM_VSC * DIM / n0 / lmd;//粘性項の計算に用いる係数
        A2 = SND * SND / n0;				//圧力の計算に用いる係数
        A3 = -DIM / n0;					//圧力勾配項の計算に用いる係数
        for (int i = 0; i < 2; i++)//粒子の種類数(0:fld,1:wll)分初期化
        {
            Dns.Add(0.0);
        }
        for (int i = 0; i < 2; i++)
        {
            invDns.Add(0.0);
        }

        Dns[FLD] = DNS_FLD;
        Dns[WLL] = DNS_WLL;
        invDns[FLD] = 1.0 / DNS_FLD;
        invDns[WLL] = 1.0 / DNS_WLL;
        rlim = PCL_DST * DST_LMT_RAT;//これ以上の粒子間の接近を許さない距離
        rlim2 = rlim * rlim;
        COL = 1.0 + COL_RAT;
        iLP = 0;            //反復数
        iF = 0;         //ファイル番号
        TIM = 0.0;		//時刻


    }

    void AlcBkt()
    {
        r = PCL_DST * 2.1;      //影響半径
        r2 = r * r;
        DB = r * (1.0 + CRT_NUM);   //バケット1辺の長さ
        DB2 = DB * DB;
        DBinv = 1.0 / DB;
        nBx = (int)((MAX_X - MIN_X) * DBinv) + 3;//解析領域内のx方向のバケット数
        nBy = (int)((MAX_Y - MIN_Y) * DBinv) + 3;//解析領域内のy方向のバケット数
        nBz = (int)((MAX_Z - MIN_Z) * DBinv) + 3;//解析領域内のz方向のバケット数
        nBxy = nBx * nBy;
        nBxyz = nBx * nBy * nBz;        //解析領域内のバケット数

        Debug.Log("nBxyz" + nBxyz);

        for (int i = 0; i < nBxyz; i++)//
        {
            bfst.Add(-1);
            blst.Add(-1);

        }
        for (int i = 0; i < nP; i++)
        {
            nxt.Add(-1);
        }
    }

    void MkBkt()//バケット格納
    {
    }
    void Vsctrm()//粘性項と重力項から仮加速度計算
    {
    }
    void UpPcl1()//仮加速度から仮速度、仮位置計算
    {
        for (int i = 0; i < nP; i++)
        {
            if (Typ[i] == 0)
            {
                Vector3 v = new Vector3(0.001f, 0.00f, 0.0f);
                Pos[i] = Pos[i] + v;
                list_particle_[i].transform.position = Pos[i];
            }
        }

    }
    void ChkCol()//剛体衝突計算
    {
    }
    void MkPrs()//仮圧力計算
    {
    }
    void PrsGrdTrm()//圧力勾配項、修正加速度計算
    {
    }
    void UpPcl2()//修正加速度から時刻更新後の位置と速度計算
    {
    }

    void ClcEMPS()
    {
        while (true)
        {
            if (iLP % 100==0)//進行表示
            {
                Debug.Log("time" + TIM + " step" + iLP);
            }
            if (TIM > FIN_TIM) break;//終了判定

            MkBkt();//バケット格納
            Vsctrm();//粘性項と重力項から仮加速度計算
            UpPcl1();//仮加速度から仮速度、仮位置計算
            ChkCol();//剛体衝突計算
            MkPrs();//仮圧力計算
            PrsGrdTrm();//圧力勾配項、修正加速度計算
            UpPcl2();////修正加速度から時刻更新後の位置と速度計算
            MkPrs();//仮圧力計算
            for (int i = 0; i < nP; i++)
            {
                pav[i] = Prs[i];
            }
            iLP++;
            TIM += DT;

         }
    }





    public GameObject Sphere;
    GameObject obj;
    // Use this for initialization
    void Start () {
        RdDat();
        
       
        for (int i = 0; i < Pos.Count(); i++)
        {
            var obj = GameObject.Instantiate(Sphere,Pos[i] , Quaternion.identity) as GameObject;
            list_particle_.Add(obj);
        }
        Debug.Log(Pos.Count());
        Debug.Log(Typ.Count());
        SetPara();
        AlcBkt();
        
    }
	
	// Update is called once per frame
	void Update () {
        //Destroy(obj);
        //UpPcl();
        //ClcEMPS();

            if (iLP % 100 == 0)//進行表示
            {
                Debug.Log("time" + TIM + " step" + iLP);
            }
            //if (TIM > FIN_TIM) break;//終了判定

            MkBkt();//バケット格納
            Vsctrm();//粘性項と重力項から仮加速度計算
            UpPcl1();//仮加速度から仮速度、仮位置計算
            ChkCol();//剛体衝突計算
            MkPrs();//仮圧力計算
            PrsGrdTrm();//圧力勾配項、修正加速度計算
            UpPcl2();////修正加速度から時刻更新後の位置と速度計算
            MkPrs();//仮圧力計算
            for (int i = 0; i < nP; i++)
            {
                pav[i] = Prs[i];
            }
            iLP++;
            TIM += DT;

        


    }




}
