using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using PacBio.Data;
using Bio;
namespace VariantCaller
{
    public class CCSRead
    {
        #region STATICS

        /// <summary>
        /// This is to avoid storing tons of unique strings with the same silly name.
        /// </summary>
        public static readonly string GENERIC_READ_NAME = "CCSRead";

        private static QualityExperiment parentExperiment = null;

        public static QualityExperiment ParentExperiment {
            get {
                return parentExperiment;
            }
            set {
                if (parentExperiment != null) {
                    throw new BioinformaticsException (
                        @"Don't set this more than once! If it is set twice.  This field exists so CCS
                          reads can get data from their parent experiment without each having a field.  If
                          there are multiple experiments in existance, please rearrange the model accordingly.");                     
                }
                parentExperiment = value;
            }
        }
        #endregion

        public readonly int ZMWnumber;
        public readonly string Movie;
        public List<CCSSubRead> SubReads;
        public Sequence Seq;
        public Zmw ZMW;

        
        /// <summary>
        /// The DNA this read is derived from.
        /// </summary>
        public Reference AssignedReference;
        /// <summary>
        /// Create a CCS read from the given sequence, usually with an ID like:
        /// "m141008_060349_42194_c100704972550000001823137703241586_s1_p0/43/ccs"
        /// </summary>
        /// <param name="read">Read.</param>
        public CCSRead(Sequence read)
        {
            //m141008_060349_42194_c100704972550000001823137703241586_s1_p0/43/ccs
            string[] sp = new string[3];
            Bio.Util.FastStringUtils.Split(read.ID,'/',sp);
            Movie = String.Intern (sp [0]);
            read.ID = GENERIC_READ_NAME;
            ZMWnumber = Convert.ToInt32(sp[1]);
            Seq = read;
        } 
    }
}
