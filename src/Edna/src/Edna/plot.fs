#light
module PacBio.Analysis.Plot


open PacBio.Common.LIMS
open PacBio.Common.IO
open PacBio.Common.Numeric
open PacBio.Common.Diagnostics

open PacBio.Analysis.LegacyDataWrapper
open Cairo


let drawPdf(run : TblRun) =

    let (width, height) = 72.*12.*20., 72.*38.
    let pdf = new PdfSurface("test.pdf", width, height)


    let ctx = new Context(pdf)
    
    //ctx.SelectFontFace("Sans", FontSlant.Normal, FontWeight.Normal)
    ctx.SetFontSize(30.)


    let sc = Seq.hd run.TblScan    
    let scan = new Scan(sc, None)    
    let trace = scan.Traces |> Seq.filter (fun tr -> tr.X = 19 && tr.Y = 38) |> Seq.hd

    let trd = trace.Trace |> Array.mapi (fun ch a -> a |> Array.map (fun i -> float(i) - trace.Baseline.[ch] ))

    
    let ymin = Stats.Percentile(trd |> Seq.concat, 0.1) * 1.3
    let ymax = Stats.Percentile(trd |> Seq.concat, 99.9) * 1.3
    let nf = trd.[0].Length
    
    let drawSeg f1 f2 yoff myheight =
    
        let t1 = float(f1) / 100.
        let t2 = float(f2) / 100.
        let t = t2-t1    

        let sx x = x / t * width
        let sy y = -yoff + height - (y - ymin) / (ymax-ymin) * myheight
        
        let dx = ref 2.1
        let dy = ref 2.1
        ctx.DeviceToUserDistance(dx, dy)
        ctx.LineWidth <- !dx
        ctx.LineJoin <- LineJoin.Round
        
        ctx.SetFontSize(30.)
        
        let textCenter x y str =
            let ex = ctx.TextExtents(str)
            let xp = (sx x) - (ex.Width/2. + ex.XBearing)
            let yp = (sy y) - (ex.Height/2. + ex.YBearing) 
            ctx.MoveTo(xp, yp)
            ctx.ShowText(str)   
        
        let setCol(r,g,b,a) = ctx.SetSourceRGBA(r,g,b,a)

        let cf = [(5., 113., 176., 200.); (26., 150., 65., 200.); (235., 110. , 0., 180.);
                    (202., 0., 32., 180.)] |> List.map (fun (a,b,c,d) -> (a/256.,b/256.,c/265., d/256.))

        let noa (a,b,c,d) = (a,b,c,1.)

        for i = 0 to 3 do
            setCol(cf.[i])
            ctx.MoveTo(0., sy (trd.[i].[0]))
            for j = f1+1 to f2-1 do
                let x = float(j-f1) / 100.
                ctx.LineTo(sx x, sy (trd.[i].[j]))
            ctx.Stroke()
            ()
            
            
        let drawPulse (p : Pulse) =
           setCol(noa cf.[p.Channel])
           textCenter (p.t1 + p.dt / 2. - t1) (p.pkmax + 200.) (sprintf "%c" p.Base)
           ctx.Stroke()
            
        trace.Pulses |> Seq.filter (fun p -> p.t1 > t1 && p.t1 < t2) |> Seq.iter drawPulse
        ()


    let n = 4
    {0..n-1} |> Seq.iter (fun i -> drawSeg (i*nf/n) ((i+1)*nf/n) (float(n-i-1)*height/float(n)) (height/float(n)))
    


    pdf.Finish() |> ignore 
    pdf.Destroy()
    ()
            
    
        
    
    
    
    
    
    
    
    