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
            int a = 0;// счетчик количества строк
            HashSet<string> samp=new HashSet<string>(); //Множество нужных людей
            List<int> samples = new List<int>();//номера столбцов vcf файла соответствующих нужным людям
            string[] splitline;
            Dictionary<string,int[]> N0= new Dictionary<string, int[]>();
            Dictionary<string, string> alleles = new Dictionary<string, string>();
            string[] readText = File.ReadAllLines("samples1000g.txt");//читаем файл с ID нужных людей
            Dictionary<char, int> stint = new Dictionary<char, int> { { '0', 0 }, { '1', 1 }, { '2', 2 }, { '3', 3 } }; //словарь для легкой конвертации цифр в числа
            foreach (string line in readText)
            {
                samp.Add(line.Split('\t')[0]);
            }
            for(int ii=0; ii < 22; ii++)//открываем файлы в цикле
            {
                samples.Clear();
                Console.WriteLine(ii);
                readText = File.ReadAllLines("Files/ALL.chr" + (ii + 1) + ".omni_2123_samples_b37_SHAPEIT.20120103.snps.chip_based.haplotypes.vcf");
                foreach (string line in readText)
                {
                    a++;
                    
                   splitline = line.Split('\t');
                    if (a == 40000)//каждые 40000 строк выводим в консоль чтобы понимать этап процесса
                    {
                        Console.WriteLine(a);
                        a = 0;
                    }
                    if (splitline.Count()>0 && splitline[0]== "#CHROM")// строка с наименованиями столбцов
                    {
                        for(int i=0; i< splitline.Count(); i++)//добавляем номера столбцов нужных людей в список samples
                        {
                            if (samp.Contains(splitline[i]))
                            {
                                samples.Add(i);
                            }
                        }
                    }
                    
                    if (splitline.Count() > 2 && splitline[1][0] <='9' && splitline[1][0] >= '0')//проверка на то, что строка подходящая
                    {
                        if (N0.Keys.Contains(splitline[0]+'_'+splitline[1])) { Console.WriteLine(splitline[0]+'_'+splitline[1]); Console.ReadLine(); } //снипы записываем как "номер хромосомы"_"позиция в хромосоме"
                        N0[splitline[0]+'_'+splitline[1]] =new int[5];
                        alleles[splitline[0]+'_'+splitline[1]] = splitline[3] + ' ' + splitline[4] + ' ' + splitline[5];//буквы нуклеотидов
                        
                        foreach(int i in samples)
                        {
                            N0[splitline[0]+'_'+splitline[1]][0] += 1;//+1 к общему числу
                            int idx = stint[splitline[i][0]];//нуклеотид на первой хромосоме
                            if (idx == 3) { Console.WriteLine(333); }
                            N0[splitline[0]+'_'+splitline[1]][idx + 1] += 1;//+1 к  числу соответствующих нуклеотидов
                            idx = stint[splitline[i][2]];//на второй хромосоме
                            N0[splitline[0]+'_'+splitline[1]][idx + 1] += 1;
                        }
                    }
                }
            }
            using (StreamWriter sw = new StreamWriter("data.txt", false, System.Text.Encoding.Default))//записываем результат
            {
                foreach(string key in N0.Keys)
                {
                    sw.WriteLine(key + ' ' + (N0[key][0]) + ' ' +(N0[key][1]) + ' ' + (N0[key][2]) + ' ' + (N0[key][3]) + ' '+ N0[key][4] + ' ' + alleles[key]);
                }
            }
        }
    }
}
