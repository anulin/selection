using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using System.IO;

namespace _1000genomes
{
    class Program
    {
        static void Main(string[] args)
        {
            Dictionary<string, int[]> N0 = new Dictionary<string, int[]>();//одна популяция
            Dictionary<string, int[]> N1 = new Dictionary<string, int[]>();// вторая популяция
            Dictionary<string, string[]> alleles1 = new Dictionary<string, string[]>();
            Dictionary<string, string[]> alleles = new Dictionary<string, string[]>();
            
            string[] readText =File.ReadAllLines("data.txt");//читаем файл первой популяции
            foreach (string line in readText)
            {
                string[] s = line.Split(' ');
                N0[s[0]] = new int[] { int.Parse(s[1])*2, int.Parse(s[2]) , int.Parse(s[3]) , int.Parse(s[4]) };//записываем частоты снипов первой популяции
                alleles[s[0]] = new string[] { s[6], s[7], s[8] };//буквы снипов
            }
            readText = File.ReadAllLines("dataC.txt");//читаем файл второй популяции
            foreach(string line in readText)
            {
                string[] s = line.Split(' ');
                N1[s[0]] = new int[] { int.Parse(s[1]) * 2, int.Parse(s[2]), int.Parse(s[3]), int.Parse(s[4]) };//записываем частоты снипов второй популяции

            }
            double n0, n1, p0, p1;//n0, n1 - общие количества образцов первой и второй популяции р0, р1 - соответствующие частоты снипов
            List<string> text = new List<string>();//список строк файла с результатами

            foreach(string key in N1.Keys)//идем по снипам
            {
                if (N1[key][0] != 0)
                {
                    n0 = N0[key][0];
                    p0 = N0[key][1] / n0;
                    n1 = N1[key][0];
                    p1 = N1[key][1] / n1;
                    if (p1 - p0 != 0)//вычисляем метрику по формуле и записываем в список результатов имя снипа, метрику, буквы
                    {
                        text.Add(key + ' ' + ((p0 - p1) * (p0 - p1) / (
                                p0 * (1 - p0) / n0 + p0 * p0 - 2 * p0 * p1 + p1 * (1 - p1) / n1 + p1 * p1 - (p1 - p0) * (p1 - p0))).ToString() + ' ' + alleles[key][0] + ' ' + alleles[key][1] + ' ' + alleles[key][2]);
                    }
                    else
                    {
                        text.Add(key + ' ' + '0' + ' ' + alleles[key][0] + ' ' + alleles[key][1] + ' ' + alleles[key][2]);
                    }
                }
            }
            File.WriteAllLines("scoresvcf.txt", text);//записываем результат в файл
            using (StreamWriter sw = new StreamWriter("tstscores.txt", false, System.Text.Encoding.Default))
            {
                sw.Write(text);
            }
            
        }
    }
}
