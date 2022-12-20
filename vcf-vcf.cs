using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using System.IO;
using System.Globalization;

namespace _1000genomes
{
    class Program
    {
        static void Main(string[] args)
        {
            double asd = 0;
            Dictionary<string, int[]> N0 = new Dictionary<string, int[]>();//одна популяция
            Dictionary<string, int[]> N1 = new Dictionary<string, int[]>();// вторая популяция
            Dictionary<string, string[]> alleles1 = new Dictionary<string, string[]>();
            Dictionary<string, string[]> alleles = new Dictionary<string, string[]>();
            Console.WriteLine(args[0]);
            int i1, i2;
            
            var maf = 0.05;
            float ind = 1;
            i1 = int.Parse(args[0]);
            i2 = int.Parse(args[1]);
            
            if (Array.IndexOf(args, "--help") != -1)
            {
                Console.Write(@"dotnet scores.dll snpcountsfile1 snpcountsfile2 samplesize1 samplesize2 [--optional flags]
--help - outputs help
--maf - minor allele frequency filter (default 0.05)
--mi - minimal fraction of individuals available for snp (default 1)");
                return;
            }
            
            int id = Array.IndexOf(args, "--maf");
            if( id!=-1){
                maf=float.Parse(args[id + 1]);
            }
            id=Array.IndexOf(args, "--mi");
            if (id != -1)
            {
                ind = float.Parse(args[id + 1]);
            }
            int muliplier = 2;
            //readText = File.ReadAllLines("syn.txt");
            System.Threading.Thread.CurrentThread.CurrentCulture = new CultureInfo("en-US");
            string[] readText =File.ReadAllLines("data.txt");//читаем файл первой популяции
            Console.WriteLine(readText[148909]);
            foreach (string line in readText)
            {

                string[] s = line.Split('\t');
                if (s[0][0] == 'Y') muliplier = 1; // Y гаплоидная
                else muliplier = 2;
                N0[s[0]] = new int[] { int.Parse(s[1])* 1, int.Parse(s[2]) , int.Parse(s[3]) , int.Parse(s[4]), int.Parse(s[5]) };//записываем частоты снипов первой популяции
                alleles[s[0]] = new string[] { s[6], s[7], s[8] };//
                
            }
            readText = File.ReadAllLines("dataC.txt");//читаем файл второй популяции
            foreach(string line in readText)
            {
                string[] s = line.Split('\t');
                if (s[0][0] == 'Y') muliplier = 1; // Y гаплоидная
                else muliplier = 2;
                N1[s[0]] = new int[] { int.Parse(s[1])* 1, int.Parse(s[2]), int.Parse(s[3]), int.Parse(s[4]), int.Parse(s[5]) };//записываем частоты снипов второй популяции

            }
            double n0, n1, p0, p1;//n0, n1 - общие количества образцов первой и второй популяции р0, р1 - соответствующие частоты снипов
            List<string> text = new List<string>();
            foreach (string key in N1.Keys)//идем по снипам
            {
                if (!N0.ContainsKey(key)) continue;
                if (N1[key][0] != 0)
                {
                    double max = 0;
                    int j = 0;
                    n0 = N0[key][0];
                    n1 = N1[key][0];
                    if (n0 <= ind * i1 || n1  <= ind * i2) continue;
                    for (int i = 1; i < 5; i++)//выбираем, которая частота нам больше подходит
                    {
                        p0 = N0[key][i] / n0;
                        p1 = N1[key][i] / n1;
                        if ((p0 - p1) * (p0 - p1) / (p0 * (1 - p0) / n0 + p1 * (1 - p1) / n1) > max)
                        {
                            max = (p0 - p1) * (p0 - p1) / (p0 * (1 - p0) / n0 + p0 * p0 - 2 * p0 * p1 + p1 * (1 - p1) / n1 + p1 * p1 - (p1 - p0) * (p1 - p0)); //(p1 - p0)^2==p0 * p0 - 2 * p0 *...
                            j = i;
                        }
                    }

                    p0 = N0[key][j] / n0;
                    p1 = N1[key][j] / n1;
                    if ((p0 <maf && p1 < maf) || ( p0 > 1- maf && p1 > 1- maf)) { continue; }//maf 0.05
                    //if(p0==p1 && (p0 == 1||p0==0)) { continue; }//отброс нулей
                    if (p1 - p0 != 0)//вычисляем метрику по формуле и записываем в список результатов имя снипа, метрику, буквы
                    {
                        text.Add(key + '\t' + max.ToString() + '\t' + alleles[key][0] + ' ' + alleles[key][1] + ' ' + alleles[key][2]);
                        if (!Double.IsPositiveInfinity(max))
                            asd += max;
                        else Console.Write("sadasdsadasdasdsada");
                    }
                    else
                    {
                        text.Add(key + '\t' + '0' + '\t' + alleles[key][0] + ' ' + alleles[key][1] + ' ' + alleles[key][2]);
                    }
                }
                else Console.WriteLine("strannoe chislo ludey");
            }
            Console.WriteLine(asd);

            File.WriteAllLines($"scoresvcf{ind}.txt", text);//записываем результат в файл
            using (StreamWriter sw = new StreamWriter("tstscores.txt", false, System.Text.Encoding.Default))
            {
                sw.Write(text);
            }
            
            Console.Read();
        }
    }
}
