#light

module PacBio.Analysis.CSV

open System.IO
open System

open Microsoft.FSharp.Control
open System.Threading
open System.Collections.Generic
open System

type ITableWriter =
    inherit IDisposable
    abstract AddIndex : (string * #seq<string>) -> unit
    abstract AddCols : (string * obj) seq -> unit
    abstract Row : unit -> unit

type IDatabase =
    inherit IDisposable
    abstract CreateTable : string -> ITableWriter
    abstract NewTransaction : unit -> unit

// 3 sig figs is enough for our purposes
let fmt i = 
    match (box i) with
    | :? float as x   -> x.ToString("F4")
    | :? float32 as x -> x.ToString("F4")
    | _               -> i.ToString().Replace(",",";")

type CSVWriter(filename) =
    let f = new FileStream(filename, FileMode.Create)
    let w = new StreamWriter(f)
    let csvLine vals = 
        w.WriteLine(String.concat "," (vals |> Seq.map fmt))
    
    let colNames = ref [| |] 
    let curRow = ref [| |]
    let columnIndicies = ref Map.empty
    let mutable firstRow = true
    
    let addCol name data =
        let n = (!curRow).Length 
        columnIndicies := (!columnIndicies).Add(name, n)
        curRow := Seq.append (!curRow) (Seq.singleton data) |> Seq.toArray
        colNames := Seq.append (!colNames) (Seq.singleton name) |> Seq.toArray
    
    
    interface ITableWriter with        
        member x.AddCols(cols : (string * obj) seq) =
            let updateCols (col, data) =
                match (!columnIndicies).TryFind(col) with
                    | Some(i) -> (!curRow).[i] <- data
                    | None -> addCol col data
            cols |> Seq.iter updateCols
            
        member x.Row() =
            if firstRow then
                csvLine (Seq.cast !colNames)
                firstRow <- false
                
            csvLine !curRow
            
        member x.AddIndex(_) = ()
            
            
    interface IDisposable with
        member x.Dispose() =
            w.Close()
            f.Close()
            w.Dispose()
            f.Dispose()       
        

let dumpCSV filename (headers: #seq<string>) (data: #seq<#seq<'a>>) =
    use f = new FileStream(filename, FileMode.Create)
    use w = new StreamWriter(f)

    let cleanString x = x.ToString().Replace(",",";")
    // Write one line to the csv file
    let csvLine vals = 
        w.WriteLine(String.concat "," (vals |> Seq.map cleanString))
        
    csvLine (Seq.cast headers)
    data |> Seq.iter (fun dataRow -> csvLine (dataRow |> Seq.cast))


let tm = [
    typeof<string>, "TEXT";
    typeof<char>, "TEXT";
    typeof<System.Double>, "NUMERIC";     
    typeof<float32>, "NUMERIC";
    typeof<int>, "NUMERIC";
    typeof<uint32>, "NUMERIC";
    typeof<uint16>, "NUMERIC"]
    
let typeMap = tm |> Seq.fold (fun (m : Map<string,string>) (typ, str) -> m.Add(typ.FullName, str)) Map.empty 
let luType (typ : Type) = if typeMap.ContainsKey(typ.FullName) then typeMap.[typ.FullName] else "TEXT"

#if FALSE

type SQLWriter(tableName, conn : SqliteConnection, trans : SqliteTransaction ref) =       
    let colNames = ref [||]
    let colTypes = ref [||]
    let curRow = ref [||]
    let pars : SqliteParameter[] ref = ref [||]
    
    let idxs : (string * seq<string>) list ref = ref []
    
    let columnIndicies = ref Map.empty
    let mutable firstRow = true
    let lastTransaction : SqliteTransaction ref = ref null
    let insertCommand : SqliteCommand ref = ref null
    
    let makeTable () =
        let cmd = conn.CreateCommand()
        cmd.Transaction <- !trans
        
        let columns = Array.zip !colNames !colTypes        
        let col (name, typ : Type) = sprintf "\"%s\" %s" name (luType typ)
        let cols = columns |> Array.map col |> (fun c -> String.Join(",", c))
        cmd.CommandText <- sprintf "CREATE TABLE \"%s\" (%s)" tableName cols
        cmd.ExecuteNonQuery()
        
    let makeIndex (name, cols) = 
        let cmd = conn.CreateCommand()
        
        let sql = sprintf "CREATE INDEX %s ON %s (%s)" name tableName (String.Join(",", cols |> Seq.map (fun s -> "[" + s + "]") |> Seq.toArray))        
        cmd.CommandText <- sql
        cmd.ExecuteNonQuery() |> ignore
        
         
    let addCol name data =
        let n = (!curRow).Length 
        columnIndicies := (!columnIndicies).Add(name, n)
        curRow := Array.append !curRow [|data|] 
        colNames := Array.append !colNames [|name|]
        colTypes := Array.append !colTypes [|data.GetType()|]
 
    
    let prepDB() =  
        makeTable() |> ignore
        !idxs |> Seq.iter makeIndex  
        let columns = Array.zip !colNames !colTypes
        pars := Array.init columns.Length (fun i -> new SqliteParameter(sprintf "@p%d" i))      
        
        
    let getCmd() = 
        if (!lastTransaction) <> (!trans) then
            lastTransaction := !trans            
            // Set up an insert command 
            let columns = Array.zip !colNames !colTypes
            let insCmd = conn.CreateCommand()
            insCmd.Transaction <- !trans
            let parNames = String.Join(",", Array.init columns.Length (sprintf "@p%d"))
            insCmd.CommandText <- sprintf "INSERT INTO %s VALUES (%s)" tableName parNames
            
            // Save the SQLite parameter for later injections                   
            insCmd.Parameters.AddRange(!pars)
            insertCommand := insCmd
            
        !insertCommand            
    
    interface ITableWriter with
        
        member x.AddCols(cols : (string * obj) seq) =
            let updateCols (col, data) =
                match (!columnIndicies).TryFind(col) with
                    | Some(i) -> (!curRow).[i] <- data
                    | None -> addCol col data
            cols |> Seq.iter updateCols
     
                  
        member x.Row() =
            if firstRow then
                prepDB()
                firstRow <- false
                
            !curRow |> Array.iteri (fun i v -> (!pars).[i].Value <- v)
            getCmd().ExecuteNonQuery() |> ignore
            
        member x.AddIndex((name, cols)) = idxs := (name, (cols :> string seq)) :: !idxs
            
            
            
    interface IDisposable with
        member x.Dispose() = ()

    
    
type SQLDb(filename) =
    let cs = sprintf "Data Source=\"%s\"" filename
    let conn = new SqliteConnection(cs)    
    let trans : SqliteTransaction ref = ref null //conn.BeginTransaction()

    let prepDB() =
        // Open the connection and a transaction 
        printfn "Connecting to: %s" cs
        conn.Open()
        trans := conn.BeginTransaction()
    do 
        prepDB()
        
        
    interface IDatabase with
        member x.CreateTable(name) = new SQLWriter(name, conn, trans) :> ITableWriter
        
        member x.NewTransaction() =
            if !trans <> null && (!trans).Connection <> null then 
                (!trans).Commit()
                trans := conn.BeginTransaction()
        
        
    interface IDisposable with
        member x.Dispose() =
            if !trans <> null && (!trans).Connection <> null then (!trans).Commit()
            conn.Close()
            conn.Dispose()        

            
let getSQLWriter filename tableName =
    if File.Exists(filename) then
        File.Delete(filename)
        
    let db = new SQLDb(filename) :> IDatabase
    db.CreateTable tableName
    
    
    
let testSql(f) =
    use w = getSQLWriter f "testtable"
    let r = new Random()
    
    let cols = ["a"; "b"; "c"; "d"]
    let rgen() = Array.init 4 (fun i -> (cols.[i], (r.NextDouble() :> obj)))
    
    Seq.init 1000 id |> Seq.iter (fun _ -> w.AddCols(rgen()); w.Row())
            


type SQLWriter(tableName, conn : SQLiteConnection, trans : SQLiteTransaction ref) =       
    let colNames = ref [||]
    let colTypes = ref [||]
    let curRow = ref [||]
    let pars : SQLiteParameter[] ref = ref [||]
    
    let idxs : (string * seq<string>) list ref = ref []
    
    let columnIndicies = ref Map.empty
    let mutable firstRow = true
    let lastTransaction : SQLiteTransaction ref = ref null
    let insertCommand : SQLiteCommand ref = ref null
    
    let makeTable () =
        let cmd = conn.CreateCommand()
        cmd.Transaction <- !trans
        
        let columns = Array.zip !colNames !colTypes        
        let col (name, typ : Type) = sprintf "\"%s\" %s" name (luType typ)
        let cols = columns |> Array.map col |> (fun c -> String.Join(",", c))
        cmd.CommandText <- sprintf "CREATE TABLE \"%s\" (%s)" tableName cols
        cmd.ExecuteNonQuery()
        
    let makeIndex (name, cols) = 
        let cmd = conn.CreateCommand()
        
        let sql = sprintf "CREATE INDEX %s ON %s (%s)" name tableName (String.Join(",", cols |> Seq.map (fun s -> "[" + s + "]") |> Seq.toArray))        
        cmd.CommandText <- sql
        cmd.ExecuteNonQuery() |> ignore
        
         
    let addCol name data =
        let n = (!curRow).Length 
        columnIndicies := (!columnIndicies).Add(name, n)
        curRow := Array.append !curRow [|data|] 
        colNames := Array.append !colNames [|name|]
        colTypes := Array.append !colTypes [|data.GetType()|]
 
    
    let prepDB() =  
        makeTable() |> ignore
        !idxs |> Seq.iter makeIndex  
        let columns = Array.zip !colNames !colTypes
        pars := Array.init columns.Length (fun i -> new SQLiteParameter(sprintf "@p%d" i))      
        
        
    let getCmd() = 
        if (!lastTransaction) <> (!trans) then
            lastTransaction := !trans            
            // Set up an insert command 
            let columns = Array.zip !colNames !colTypes
            let insCmd = conn.CreateCommand()
            insCmd.Transaction <- !trans
            let parNames = String.Join(",", Array.init columns.Length (sprintf "@p%d"))
            insCmd.CommandText <- sprintf "INSERT INTO %s VALUES (%s)" tableName parNames
            
            // Save the SQLite parameter for later injections                   
            insCmd.Parameters.AddRange(!pars)
            insertCommand := insCmd
            
        !insertCommand            
    
    interface ITableWriter with
        
        member x.AddCols(cols : (string * obj) seq) =
            let updateCols (col, data) =
                match (!columnIndicies).TryFind(col) with
                    | Some(i) -> (!curRow).[i] <- data
                    | None -> addCol col data
            cols |> Seq.iter updateCols
     
                  
        member x.Row() =
            if firstRow then
                prepDB()
                firstRow <- false
                
            !curRow |> Array.iteri (fun i v -> (!pars).[i].Value <- v)
            getCmd().ExecuteNonQuery() |> ignore
            
        member x.AddIndex((name, cols)) = idxs := (name, (cols :> string seq)) :: !idxs
            
            
            
    interface IDisposable with
        member x.Dispose() = ()





    
    
type SQLDb(filename) =
    let conn = new SQLiteConnection(sprintf "Data Source=\"%s\"" filename)    
    let trans : SQLiteTransaction ref = ref null //conn.BeginTransaction()

    let prepDB() =
        // Open the connection and a transaction 
        conn.Open()
        trans := conn.BeginTransaction()
    do 
        prepDB()
        
        
    interface IDatabase with
        member x.CreateTable(name) = new SQLWriter(name, conn, trans) :> ITableWriter
        
        member x.NewTransaction() =
            if !trans <> null && (!trans).Connection <> null then 
                (!trans).Commit()
                trans := conn.BeginTransaction()
        
        
    interface IDisposable with
        member x.Dispose() =
            if !trans <> null && (!trans).Connection <> null then (!trans).Commit()
            conn.Close()
            conn.Dispose()        

            
let getSQLWriter filename tableName =
    if File.Exists(filename) then
        File.Delete(filename)
        
    let db = new SQLDb(filename) :> IDatabase
    db.CreateTable tableName
    
    
    
let testSql(f) =
    use w = getSQLWriter f "testtable"
    let r = new Random()
    
    let cols = ["a"; "b"; "c"; "d"]
    let rgen() = Array.init 4 (fun i -> (cols.[i], (r.NextDouble() :> obj)))
    
    Seq.init 1000 id |> Seq.iter (fun _ -> w.AddCols(rgen()); w.Row())

#endif
