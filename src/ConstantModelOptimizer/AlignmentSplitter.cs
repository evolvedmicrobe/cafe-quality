using System;
using System.Collections.Generic;
using System.Linq;
using Bio;
using Bio.Algorithms.Alignment;
using VariantCaller;

namespace ConstantModelOptimizer
{
    // Template then read position.
    using Positions = KeyValuePair<int, int>;

    public class AlignmentSplitter
    {
        const byte GAP_CHAR = (byte)'-';
        /// <summary>
        /// This method takes a read and template and attempts to divide it into smaller fragments for dynamic programming.
        /// It does this by taking a MUMMER alignment of the two sequences and then identifying 'anchor points' as places where 
        /// 5 basepairs continuously match.  It then creates sub sequences between anchor points separated by approximately 50 template
        /// bases.
        /// </summary>
        /// <returns>The sequence.</returns>
        /// <param name="read">Read.</param>
        /// <param name="template">Template.</param>
        public static List<ReadTemplatePair> SplitSequence(string read, string template)
        {
            try {
            var zzz = read;

            // Get the template and the read sequence
            var tpl = new Sequence(NoGapDnaAlphabet.Instance, template, false);
            var rd = new Sequence(NoGapDnaAlphabet.Instance, read, false);
            var nucmer= new NucmerQueryable(tpl, 7);

            // Align them
            var toR = nucmer.GetAlignments(rd);
            toR.Sort ((x, y) => -x.Score.CompareTo (y.Score));
//            Console.WriteLine (template);
//            Console.WriteLine (read);
//            Console.WriteLine ("Start Length: "+read.Length);
            // Simple heuristic, we will define a plausible break as any point in the alignment
            // where it seems like the center basepair and the two on either side are connected
            // Let's make a list of all such positions
            List<Positions> matches = new List<Positions> ();
            foreach (var aln in toR) {
                var t = ((Sequence)aln.FirstSequence).GetInternalArray();
                var r = ((Sequence)aln.SecondSequence).GetInternalArray();
                var currentMatches = 0;
                var temp_pos = 0;
                var read_pos = 0;
                //Console.WriteLine (aln);
                for (int i = 0; i < t.Length; i++) {
                    var tbp = t[i];
                    var rbp = r [i];
                    var gt = tbp != GAP_CHAR;
                    var gr = tbp != GAP_CHAR;
                    var isMatch = tbp==rbp && gt;
                    if (isMatch) {
                        currentMatches++;
                    } else {
                        currentMatches = 0;
                    }
                    if (currentMatches >= 5) {
                        var spot_t = temp_pos - 2 + aln.FirstSequenceStart.Value;
                        var spot_r = read_pos - 2 + aln.SecondSequenceStart.Value;
                        //Console.WriteLine (spot_t+"\t"+ spot_r);
                        matches.Add (new Positions ((int)spot_t, (int)spot_r));
                    }
                    if (gt) {
                        temp_pos++;
                    }
                    if (gr) {
                        read_pos++;
                    }
                }
            }

            // Now to sort and filter the list
            matches.Sort ((x, y) => x.Key.CompareTo (y.Key));
            if (matches.First ().Key > matches.Last ().Key) {
                throw new Exception ();
            }

            // Let's curated it down to items within ~50 bp of each other
            // And simulataneously make sure we monotonically increase
            List<Positions> newPositions = new List<Positions> (20);
            var lastTemplatePos = 0;
            var idealDistance = 50;
            var lastTemp = -1;
            var lastRead = -1;
            foreach (var pos in matches) {
                if (pos.Key < lastTemp || pos.Value < lastRead) {
                    // totally poorly formatted.
                    return null;
                }
                lastTemp = pos.Key;
                lastRead = pos.Value;
                if (pos.Key > (idealDistance + lastTemplatePos)) {
                    newPositions.Add (pos);
                    lastTemplatePos = pos.Key;
                }
            }

            // Now to divy up into sequences based on these anchor points
            List<ReadTemplatePair> pairs = new List<ReadTemplatePair> ();
            int lastTemplate = 0;
            lastRead = 0;
//            Console.WriteLine ("Later Length: "+read.Length);
//            Console.WriteLine ("Later Length 2: "+zzz.Length);

                foreach (var pos in newPositions.Take(newPositions.Count -1 )) {
                var tpl_sec = template.Substring(lastTemplate, pos.Key - lastTemplate + 1);
                var rd_sec = read.Substring (lastRead, pos.Value - lastRead + 1);
                lastTemplate = pos.Key + 1;
                lastRead = pos.Value + 1;
                pairs.Add(new ReadTemplatePair(rd_sec, tpl_sec));
                }

            if (lastTemplate < template.Length && lastRead < read.Length) {
                var tpl_sec = template.Substring(lastTemplate, (template.Length-1) - lastTemplate + 1);
                var rd_sec = read.Substring (lastRead, (read.Length-1) - lastRead + 1);
                pairs.Add(new ReadTemplatePair(rd_sec, tpl_sec));
            } else {
                throw new Exception ("You now need to handle the end case");
            }
            var newTpl = pairs.Aggregate ("", (x, y) => x + y.Template);
            var newRead = pairs.Aggregate ("", (x, y) => x + y.Read);
            if (newTpl != template || newRead != read) {
                throw new Exception ("Mismatch");
            }
            return pairs;
            }
            catch(Exception thrown) {
                return null;
            }
        }
    }
}

