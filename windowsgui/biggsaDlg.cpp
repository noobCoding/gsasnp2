
// biggsaDlg.cpp : implementation file
//

#include "stdafx.h"
#include "biggsa.h"
#include "biggsaDlg.h"
#include "afxdialogex.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

initparam initParam;
resparam  resParam;

CString inputFile;
CString mapFile;
CString setFile;

bool	setWrite = 0;
bool	inputType = SNP_INPUT;

struct mythreadinitparam {
	CBigGsaDlg	*thisDlg;
	initparam	*initParam;

	virtual ~mythreadinitparam() {};
};

// CAboutDlg dialog used for App About

class CAboutDlg : public CDialogEx
{
public:
	CAboutDlg();

	// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_ABOUTBOX };
#endif

protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

// Implementation
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialogEx(IDD_ABOUTBOX)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialogEx)
END_MESSAGE_MAP()

// CBigGsaDlg dialog


CBigGsaDlg::CBigGsaDlg(CWnd* pParent /*=NULL*/)
	: CDialogEx(IDD_BIGGSA_DIALOG, pParent)
	, m_rad_hgNumber(0)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDI_ICON1);
}

void CBigGsaDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	DDX_Control(pDX, IDC_LIST_Result, m_cListCtrl);
	DDX_Control(pDX, IDC_Tab_ListType, m_cTabListMode);
	DDX_Control(pDX, IDC_LOG, m_Log);
	DDX_Control(pDX, IDC_MFCB_SetInfo, m_fSetInfo);
	DDX_Control(pDX, IDC_MFCB_GeneMap, m_fGeneMap);
	DDX_Control(pDX, IDC_MFCB_INPUT, m_fInput);
	DDX_Control(pDX, IDC_Radio_hg18, m_cRadio_hg18);
	DDX_Control(pDX, IDC_Combo_Padding, m_cCombo_Padding);
	DDX_Control(pDX, IDC_Radio_GO, m_cRadio_GO);
	DDX_Control(pDX, IDC_Combo_MSigDB, m_cCombo_MSigDB);
	DDX_Control(pDX, IDC_MinGeneSize, m_cEdit_MinGeneSize);
	DDX_Control(pDX, IDC_MaxGeneSize, m_cEdit_MaxGeneSize);
	DDX_Control(pDX, IDC_btnPerform, m_btnPerform);
	DDX_Control(pDX, IDC_MFCB_OutputFile, m_EBC_OutputFile);
	DDX_Control(pDX, IDC_ToSymbol, m_cbt_ToSymbol);
	DDX_Control(pDX, IDC_COMBO2, m_Cmb_GeneCutoff);
	DDX_Control(pDX, IDC_COMBO1, m_Cmb_AdjacentGeneFile);
	DDX_Control(pDX, IDC_CoreNet, m_CoreNet);
	DDX_Control(pDX, IDC_COMBO4, m_cmb_Net_qvalue_cutoff);
	DDX_Control(pDX, IDC_COMBO5, m_cmb_Net_GeneScoreCutoff);
	DDX_Control(pDX, IDC_RADIO_SNP, m_rad_Input_SNP);
	DDX_Control(pDX, IDC_RADIO_GENE, m_rad_Input_GENE);	
	DDX_Control(pDX, IDC_MFCB_NetFile, m_cbrowser_NetFile);
}

BEGIN_MESSAGE_MAP(CBigGsaDlg, CDialogEx)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	ON_BN_CLICKED(IDC_Radio_hg18, &CBigGsaDlg::OnBnClickedRadiohg18)
	ON_BN_CLICKED(IDC_Radio_hgUser, &CBigGsaDlg::OnBnClickedRadiohguser)
	ON_BN_CLICKED(IDC_Radio_GO, &CBigGsaDlg::OnBnClickedRadioGo)
	ON_BN_CLICKED(IDC_Radio_Other, &CBigGsaDlg::OnBnClickedRadioOther)
	ON_BN_CLICKED(IDC_Radio_KEGG, &CBigGsaDlg::OnBnClickedRadioKegg)
	ON_BN_CLICKED(IDC_Radio_MSigDB, &CBigGsaDlg::OnBnClickedRadioMsigdb)
	ON_BN_CLICKED(IDC_Radio_hg38, &CBigGsaDlg::OnBnClickedRadiohg38)
	ON_BN_CLICKED(IDC_Radio_hg19, &CBigGsaDlg::OnBnClickedRadiohg19)
	ON_NOTIFY(HDN_ITEMCLICK, 0, &CBigGsaDlg::OnHdnItemclickListResult)
	ON_BN_CLICKED(IDC_btnPerform, &CBigGsaDlg::OnBnClickedbtnperform)
	ON_EN_UPDATE(IDC_MFCB_OutputFile, &CBigGsaDlg::OnEnUpdateMfcbOutputfile)
	ON_NOTIFY(LVN_ITEMCHANGED, IDC_LIST_Result, &CBigGsaDlg::OnLvnItemchangedListResult)
	ON_NOTIFY(NM_CUSTOMDRAW, IDC_LIST_Result, &CBigGsaDlg::OnNMCustomdrawListResult)
	ON_BN_CLICKED(IDC_CoreNet, &CBigGsaDlg::OnBnClickedCorenet)
	ON_BN_CLICKED(IDC_RADIO_SNP, &CBigGsaDlg::OnBnClickedRadioSnp)
	ON_BN_CLICKED(IDC_RADIO_GENE, &CBigGsaDlg::OnBnClickedRadioGene)	
END_MESSAGE_MAP()

// CBigGsaDlg message handlers

BOOL CBigGsaDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// Add "About..." menu item to system menu.

	// IDM_ABOUTBOX must be in the system command range.
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		BOOL bNameValid;
		CString strAboutMenu;
		bNameValid = strAboutMenu.LoadString(IDS_ABOUTBOX);
		ASSERT(bNameValid);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// Set the icon for this dialog.  The framework does this automatically
	//  when the application's main window is not a dialog
	HICON m_hIcon = LoadIcon(AfxGetInstanceHandle(), MAKEINTRESOURCE(IDI_ICON2));
	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, TRUE);		// Set small icon

	// Initial extended style for the list control on this dialog
	DWORD dwStyle = m_cListCtrl.GetExtendedStyle();
	dwStyle |= LVS_EX_GRIDLINES | LVS_EX_FULLROWSELECT;
	m_cListCtrl.SetExtendedStyle(dwStyle);

	// Setup the tab header
	InitTabCtrl();

	// default input
	m_fInput.SetWindowTextW(_T("data\\DIAGRAM"));
	inputFile = _T("data\\DIAGRAM");
	m_EBC_OutputFile.SetWindowTextW(_T("DIAGRAM_zscore.txt"));
	m_cEdit_MinGeneSize.SetWindowTextW(_T("10"));
	m_cEdit_MaxGeneSize.SetWindowTextW(_T("200"));
	m_cCombo_Padding.AddString(_T("±20000"));
	m_cCombo_Padding.AddString(_T("±10000"));
	m_cCombo_Padding.AddString(_T("0"));
	m_cCombo_Padding.AddString(_T("Exon only"));
	m_Cmb_GeneCutoff.SetWindowTextW(_T("0.05"));
	m_Cmb_AdjacentGeneFile.AddString(_T("ALL races"));
	m_Cmb_AdjacentGeneFile.AddString(_T("African"));
	m_Cmb_AdjacentGeneFile.AddString(_T("American"));
	m_Cmb_AdjacentGeneFile.AddString(_T("East Asian"));
	m_Cmb_AdjacentGeneFile.AddString(_T("European"));
	m_Cmb_AdjacentGeneFile.AddString(_T("South Asian"));
	m_Cmb_AdjacentGeneFile.SetWindowTextW(_T("European"));
	m_cCombo_MSigDB.AddString(_T("c1.all.v5.2.symbols.gmt"));
	m_cCombo_MSigDB.AddString(_T("c2.all.v5.2.symbols.gmt"));
	m_cCombo_MSigDB.AddString(_T("c2.cp.v5.2.symbols.gmt"));
	m_cCombo_MSigDB.AddString(_T("c5.all.v5.2.symbols.gmt"));
	m_cmb_Net_GeneScoreCutoff.AddString(_T("0.05"));
	m_cmb_Net_GeneScoreCutoff.AddString(_T("0.01"));
	m_cmb_Net_GeneScoreCutoff.AddString(_T("0.005"));
	m_cmb_Net_GeneScoreCutoff.AddString(_T("0.001"));
	m_cmb_Net_GeneScoreCutoff.SetWindowTextW(_T("0.001"));
	m_cmb_Net_qvalue_cutoff.AddString(_T("0.05"));
	m_cmb_Net_qvalue_cutoff.AddString(_T("0.10"));
	m_cmb_Net_qvalue_cutoff.AddString(_T("0.15"));
	m_cmb_Net_qvalue_cutoff.AddString(_T("0.20"));
	m_cmb_Net_qvalue_cutoff.AddString(_T("0.25"));
	m_cmb_Net_qvalue_cutoff.SetWindowTextW(_T("0.10"));
	m_rad_Input_SNP.SetCheck(1);	
	m_cbrowser_NetFile.SetWindowTextW(_T("data\\STRING_NETWORK.txt"));
	//////////////////////
	col_Index = 6;
	method_Index = 0;
	m_truncation_ratio = 0.25;
	b_Convert2Symbol = false;
	UpdateData();
	return TRUE;  // return TRUE  unless you set the focus to a control
}

void CBigGsaDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialogEx::OnSysCommand(nID, lParam);
	}
}

void CBigGsaDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // device context for painting

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// Center icon in client rectangle
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Draw the icon
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialogEx::OnPaint();
	}
}

HCURSOR CBigGsaDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}

UINT Init8(void* pParam)
{
	CBigGsaDlg* tmp = (CBigGsaDlg*)pParam;
	CString mess;

	mess.Format(_T("Loading input data..."));
	tmp->m_Log.AppendString(mess);

	chrono::time_point<chrono::system_clock> start, end;
	start = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds;

	/*tmp->test->snpPvalueLoad(initParam.inputFile);
	tmp->test->snpGeneMapGenerate(initParam.snpGeneMapFile, initParam.inputType);	
	tmp->test->setGeneLoad(initParam.setFile, initParam.convert2Symbol);*/
	tmp->test->Initializing(initParam);

	end = chrono::system_clock::now();
	elapsed_seconds = end - start;
	mess.Format(_T("Data initialized in : %fs."), elapsed_seconds.count());
	tmp->m_Log.AppendString(mess);

	//////////////////////////////////////////////////////////////////////////////////
	mess.Format(_T("Analyzing..."));

	tmp->m_Log.AppendString(mess);
	start = chrono::system_clock::now();

	tmp->test->Analyzing(initParam);	//resParam);

	end = chrono::system_clock::now();
	elapsed_seconds = end - start;
	mess.Format(_T("Analyzed in: %f s."), elapsed_seconds.count());
	tmp->m_Log.AppendString(mess);

	//////////////////////////////////////////////////////////////////////////////////
	mess.Format(_T("Network is being loaded."));
	tmp->m_Log.AppendString(mess);
	
	////
	tmp->InitListCtrlSetView8();

	// ********************************************************
	tmp->m_btnPerform.SetWindowTextW(_T("N E T W O R K  I S  B E I N G  L O A D E D"));		
	tmp->test->networkLoading(initParam.networkFile);
	//string tmpstr = tmp->test->commonNetworkLoading(); // assuming networkLoading("data//STRING_NETWORK.txt");
	// ********************************************************
	mess.Format(_T("All done."));
	tmp->m_Log.AppendString(mess);
	mess.Format(_T("					"));
	tmp->m_Log.AppendString(mess);

	//// GUI manipulation ////////////////////////////////////////////////////////////		
	tmp->m_cTabListMode.EnableWindow(1);
	tmp->m_cListCtrl.EnableWindow(1);
	tmp->m_fInput.EnableWindow(1);
	tmp->m_cEdit_MaxGeneSize.EnableWindow(1);
	tmp->m_cEdit_MinGeneSize.EnableWindow(1);
	tmp->m_btnPerform.EnableWindow(1);
	tmp->m_EBC_OutputFile.EnableWindow(1);
	tmp->m_cCombo_Padding.EnableWindow(1);
	tmp->m_cCombo_MSigDB.EnableWindow(1);
	tmp->m_Cmb_GeneCutoff.EnableWindow(1);
	tmp->m_Cmb_AdjacentGeneFile.EnableWindow(1);
	tmp->m_CoreNet.EnableWindow(1);
	tmp->m_cmb_Net_qvalue_cutoff.EnableWindow(1);
	tmp->m_cmb_Net_GeneScoreCutoff.EnableWindow(1);
	tmp->m_cbrowser_NetFile.EnableWindow(1);
	tmp->m_btnPerform.SetWindowTextW(_T("R U N"));

	//AfxEndThread((UINT)0, 1);
	return 0;
}

UINT Init9(void* pParam) // for gene
{
	CBigGsaDlg* tmp = (CBigGsaDlg*)pParam;
	CString mess;

	mess.Format(_T("Loading input data..."));
	tmp->m_Log.AppendString(mess);

	chrono::time_point<chrono::system_clock> start, end;
	start = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds;

	CT2CA pszName2(inputFile);
	string string2File2(pszName2);
	tmp->test->genePvalueLoad(string2File2);

	mess.Format(_T("Pre-processing input data..."));
	tmp->m_Log.AppendString(mess);

	CT2CA pszName(mapFile);
	string string2File(pszName);
	tmp->test->snpGeneMapGenerate(string2File, GENE_INPUT);

	CT2CA pszName1(setFile);
	string string2File1(pszName1);
	tmp->test->setGeneLoad(string2File1, tmp->b_Convert2Symbol);

	end = chrono::system_clock::now();
	elapsed_seconds = end - start;
	mess.Format(_T("Data initialized in : %fs."), elapsed_seconds.count());
	tmp->m_Log.AppendString(mess);

	//////////////////////////////////////////////////////////////////////////////////
	mess.Format(_T("Analyzing..."));

	tmp->m_Log.AppendString(mess);
	start = chrono::system_clock::now();

	tmp->test->geneScoring(GENE_INPUT); // original
	tmp->test->adjustedGeneScoring(GENE_INPUT); // original
	tmp->test->adjustedPvalueBH();

	end = chrono::system_clock::now();
	elapsed_seconds = end - start;
	mess.Format(_T("Analysed in: %f s."), elapsed_seconds.count());
	tmp->m_Log.AppendString(mess);

	//////////////////////////////////////////////////////////////////////////////////
	mess.Format(_T("Network is being loaded."));
	tmp->m_Log.AppendString(mess);
	
	tmp->InitListCtrlSetView8();
	tmp->m_cTabListMode.EnableWindow(1);
	// ********************************************************
	tmp->m_btnPerform.SetWindowTextW(_T("N E T W O R K  I S  B E I N G  L O A D E D"));

	CT2CA psznetdata(tmp->m_networkdata);
	string strnetdata(psznetdata);
	tmp->test->networkLoading(strnetdata);
	
	// ********************************************************	
	mess.Format(_T("All done."));
	tmp->m_Log.AppendString(mess);
	mess.Format(_T("					"));
	tmp->m_Log.AppendString(mess);	

	//// GUI manipulation ////////////////////////////////////////////////////////////	
	tmp->m_cListCtrl.EnableWindow(1);
	tmp->m_fInput.EnableWindow(1);
	tmp->m_cEdit_MaxGeneSize.EnableWindow(1);
	tmp->m_cEdit_MinGeneSize.EnableWindow(1);
	tmp->m_btnPerform.EnableWindow(1);
	tmp->m_EBC_OutputFile.EnableWindow(1);
	tmp->m_cCombo_Padding.EnableWindow(1);
	tmp->m_cCombo_MSigDB.EnableWindow(1);
	tmp->m_Cmb_GeneCutoff.EnableWindow(1);
	tmp->m_Cmb_AdjacentGeneFile.EnableWindow(1);
	tmp->m_CoreNet.EnableWindow(1);
	tmp->m_cmb_Net_qvalue_cutoff.EnableWindow(1);
	tmp->m_cmb_Net_GeneScoreCutoff.EnableWindow(1);
	tmp->m_cbrowser_NetFile.EnableWindow(1);
	tmp->m_btnPerform.SetWindowTextW(_T("R U N"));

	//AfxEndThread((UINT)0, 1);
	return 0;
}

void CBigGsaDlg::InsertSetItems8()
{
	// Delete the current contents
	m_cListCtrl.DeleteAllItems();

	const int szSetInput = test->setGlobal.size();
	map<double, setinfo> mySet;
	double myepsilon = 0.000000001;
	int k = 0; // k * epsilon

	for (auto& vi : test->setGlobal)
	{
		setinfo tmp;
		tmp.setname = vi;
		tmp.genecount = test->setGeneGlobal.operator[](tmp.setname).size();
		tmp.setsize = test->setGlobalActualSize.operator[](tmp.setname);
		tmp.zscore = test->zsOrigin.operator[](tmp.setname);
		tmp.adjZscore = test->zsAdjusted.operator[](tmp.setname);
		tmp.adjPval = test->pvAdjusted.operator[](tmp.setname);
		tmp.adjQval = test->fdrAdjusted.operator[](tmp.setname);

		switch (col_Index)
		{
		case 0:
		case 1:
			//mySet[tmp.genecount] = tmp;  // sorting based on  gene count			
		case 2:
			//mySet[tmp.setsize] = tmp;  // sorting based on  set size			
		case 3:
			if (mySet.find(tmp.zscore) == mySet.end())
				mySet[tmp.zscore] = tmp;  // sorting based on  original zscore		
			else
				mySet[tmp.zscore + ++k * myepsilon] = tmp;  // sorting based on  original zscore		
			break;
		case 4:
			if (mySet.find(tmp.adjZscore) == mySet.end())
				mySet[tmp.adjZscore] = tmp; // adjusted z-score;
			else
				mySet[tmp.adjZscore + ++k * myepsilon] = tmp; //adjusted z-score;
			break;
		case 5:
			if (mySet.find(tmp.adjPval) == mySet.end()) // not exist
				mySet[tmp.adjPval] = tmp; // adjusted p-value
			else
				mySet[tmp.adjPval + ++k * myepsilon] = tmp; // adjusted p-value		
		case 7:
		case 8:
		case 6:
			if (mySet.find(tmp.adjQval) == mySet.end()) // not exist
				mySet[tmp.adjQval] = tmp; // adjusted q-value
			else
				mySet[tmp.adjQval + ++k * myepsilon] = tmp; // adjusted q-value
			break;
		}
	}

	// Use the LV_ITEM structure to insert the items
	LVITEM lvi;
	int idx = 0;
	CString strItem;

	if (col_Index != 3 && col_Index != 4)
	{
		map<double, setinfo>::reverse_iterator it;
		for (it = mySet.rbegin(); it != mySet.rend(); ++it)
			//for (auto& mi : mySet)
		{
			// Insert the first item
			lvi.mask = LVIF_IMAGE | LVIF_TEXT;
			strItem.Format(_T("%S"), it->second.setname.c_str());

			lvi.iItem = idx;
			lvi.iSubItem = 0;
			lvi.pszText = (LPTSTR)(LPCTSTR)(strItem);
			lvi.iImage = idx % 8;		// There are 8 images in the image list
			m_cListCtrl.InsertItem(&lvi);

			// Set subitem 1
			strItem.Format(_T("%d"), it->second.setsize);
			lvi.iSubItem = 1;
			lvi.pszText = (LPTSTR)(LPCTSTR)(strItem);
			m_cListCtrl.SetItem(&lvi);

			// Set subitem 2
			strItem.Format(_T("%d"), it->second.genecount);
			lvi.iSubItem = 2;
			lvi.pszText = (LPTSTR)(LPCTSTR)(strItem);
			m_cListCtrl.SetItem(&lvi);

			// Set subitem 3
			strItem.Format(_T("%.8f"), it->second.zscore);
			lvi.iSubItem = 3;
			lvi.pszText = (LPTSTR)(LPCTSTR)(strItem);
			m_cListCtrl.SetItem(&lvi);

			// Set subitem 4
			strItem.Format(_T("%.8f"), it->second.adjZscore);
			lvi.iSubItem = 4;
			lvi.pszText = (LPTSTR)(LPCTSTR)(strItem);
			m_cListCtrl.SetItem(&lvi);

			// Set subitem 5
			strItem.Format(_T("%.8f"), it->second.adjPval);
			lvi.iSubItem = 5;
			lvi.pszText = (LPTSTR)(LPCTSTR)(strItem);
			m_cListCtrl.SetItem(&lvi);

			// Set subitem 6
			strItem.Format(_T("%.8f"), it->second.adjQval);
			lvi.iSubItem = 6;
			lvi.pszText = (LPTSTR)(LPCTSTR)(strItem);
			m_cListCtrl.SetItem(&lvi);

			// Set subitem 7
			strItem.Format(_T("Show"));	//, it->second.adjQval);
			lvi.iSubItem = 7;
			lvi.pszText = (LPTSTR)(LPCTSTR)(strItem);
			m_cListCtrl.SetItem(&lvi);

			// Set subitem 8
			strItem.Empty();
			CString tmp;

			k = 1; // k * epsilon
			map<double, string> geneScore; // sort gene score & present them
			for (auto& vi : test->setGeneRefined[it->second.setname]) {
				double tmpscore = test->adjGeneScore[test->geneInputIndex[vi]];
				if (geneScore.find(tmpscore) != geneScore.end()) { // exist score
					geneScore[tmpscore + ++k * myepsilon] = vi;
				}
				else {
					geneScore[tmpscore] = vi;
				}
			}

			map<double, string>::reverse_iterator ri;
			for (ri = geneScore.rbegin(); ri != geneScore.rend(); ri++)
			{
				tmp.Format(_T("%S (%f); "), ri->second.c_str(), ri->first);
				strItem += tmp;
			}

			lvi.iSubItem = 8;
			lvi.pszText = (LPTSTR)(LPCTSTR)(strItem);
			m_cListCtrl.SetItem(&lvi);
		}
	}
	else {
		map<double, setinfo>::iterator rit;
		for (rit = mySet.begin(); rit != mySet.end(); ++rit)
		{
			// Insert the first item
			lvi.mask = LVIF_IMAGE | LVIF_TEXT;
			strItem.Format(_T("%S"), rit->second.setname.c_str());

			lvi.iItem = idx;
			lvi.iSubItem = 0;
			lvi.pszText = (LPTSTR)(LPCTSTR)(strItem);
			lvi.iImage = idx % 8;		// There are 8 images in the image list
			m_cListCtrl.InsertItem(&lvi);

			// Set subitem 1
			strItem.Format(_T("%d"), rit->second.setsize);
			lvi.iSubItem = 1;
			lvi.pszText = (LPTSTR)(LPCTSTR)(strItem);
			m_cListCtrl.SetItem(&lvi);

			// Set subitem 2
			strItem.Format(_T("%d"), rit->second.genecount);
			lvi.iSubItem = 2;
			lvi.pszText = (LPTSTR)(LPCTSTR)(strItem);
			m_cListCtrl.SetItem(&lvi);

			// Set subitem 3
			strItem.Format(_T("%.8f"), rit->second.zscore);
			lvi.iSubItem = 3;
			lvi.pszText = (LPTSTR)(LPCTSTR)(strItem);
			m_cListCtrl.SetItem(&lvi);

			// Set subitem 4			
			strItem.Format(_T("%.8f"), rit->second.adjZscore);
			lvi.iSubItem = 4;
			lvi.pszText = (LPTSTR)(LPCTSTR)(strItem);
			m_cListCtrl.SetItem(&lvi);

			// Set subitem 5
			strItem.Format(_T("%.8f"), rit->second.adjPval);
			lvi.iSubItem = 5;
			lvi.pszText = (LPTSTR)(LPCTSTR)(strItem);
			m_cListCtrl.SetItem(&lvi);

			// Set subitem 6
			strItem.Format(_T("%.8f"), rit->second.adjQval);
			lvi.iSubItem = 6;
			lvi.pszText = (LPTSTR)(LPCTSTR)(strItem);
			m_cListCtrl.SetItem(&lvi);

			// Set subitem 7
			strItem.Format(_T("Show")); 	//, rit->second.adjQval);
			lvi.iSubItem = 7;
			lvi.pszText = (LPTSTR)(LPCTSTR)(strItem);
			m_cListCtrl.SetItem(&lvi);

			// Set subitem 8
			strItem.Empty();
			CString tmp;

			k = 1; // k * epsilon
			map<double, string> geneScore; // sort gene score & present them

			for (auto& vi : test->setGeneRefined[rit->second.setname]) {
				double tmpscore = test->adjGeneScore[test->geneInputIndex[vi]];
				if (geneScore.find(tmpscore) != geneScore.end()) { // exist score
					geneScore[tmpscore + ++k * myepsilon] = vi;
				}
				else {
					geneScore[tmpscore] = vi;
				}
			}


			map<double, string>::reverse_iterator ri;
			for (ri = geneScore.rbegin(); ri != geneScore.rend(); ri++)
			{
				tmp.Format(_T("%S (%f); "), ri->second.c_str(), ri->first);
				strItem += tmp;
			}

			lvi.iSubItem = 8;
			lvi.pszText = (LPTSTR)(LPCTSTR)(strItem);
			m_cListCtrl.SetItem(&lvi);
		}
	}

	if (!setWrite) {
		map<double, setinfo> set2File;
		double kset = 0; // k * epsilon

		for (auto& vi : test->setGlobal)
		{
			setinfo tmp;
			tmp.setname = vi;
			tmp.genecount = test->setGeneGlobal.operator[](tmp.setname).size();
			tmp.setsize = test->setGlobalActualSize.operator[](tmp.setname);
			tmp.zscore = test->zsOrigin.operator[](tmp.setname);
			tmp.adjZscore = test->zsAdjusted.operator[](tmp.setname);
			tmp.adjPval = test->pvAdjusted.operator[](tmp.setname);
			tmp.adjQval = test->fdrAdjusted.operator[](tmp.setname);

			if (set2File.find(tmp.adjQval) == mySet.end())
				set2File[tmp.adjQval] = tmp;  // sorting based on  adjusted zscore		
			else
				set2File[tmp.adjQval + ++kset * myepsilon] = tmp;  // sorting based on  original zscore		
		}

		// printing output	
		CString outputfile;
		m_EBC_OutputFile.GetWindowText(outputfile);
		if (outputfile.IsEmpty()) outputfile = "adjusted_scores.txt";

		ofstream fou(outputfile);
		fou << "Set\tSize\tCount\tz-score\tAdj. z-score\tAdj. p-value\tq-value\tList of Genes\n";
		map<double, setinfo>::iterator rit;
		for (rit = set2File.begin(); rit != set2File.end(); ++rit)
		{
			// Insert the first item
			fou << rit->second.setname << "\t";
			fou << rit->second.setsize << "\t";
			fou << rit->second.genecount << "\t";
			fou << rit->second.zscore << "\t";
			fou << rit->second.adjZscore << "\t";
			fou << rit->second.adjPval << "\t";
			fou << rit->second.adjQval << "\t";

			int kgene = 0; // k * epsilon
			map<double, string> geneScore; // sort gene score & present them

			for (auto& vi : test->setGeneRefined[rit->second.setname]) {
				double tmpscore = test->adjGeneScore[test->geneInputIndex[vi]];
				if (geneScore.find(tmpscore) != geneScore.end()) { // exist score
					geneScore[tmpscore + ++k * myepsilon] = vi;
				}
				else {
					geneScore[tmpscore] = vi;
				}
			}

			map<double, string>::reverse_iterator ri;
			for (ri = geneScore.rbegin(); ri != geneScore.rend(); ri++)
			{
				fou << ri->second.c_str() << "(" << ri->first << "); ";
			}
			fou << endl;
		}
		fou.close();
		setWrite = 1;
	}
}

void CBigGsaDlg::InitListCtrlSetView8()
{
	// Insert some columns
	CRect rect;
	m_cListCtrl.GetClientRect(&rect);
	m_cListCtrl.SetRedraw(FALSE);

	// Remove whatever style is there currently
	m_cListCtrl.ModifyStyle(LVS_ICON | LVS_LIST | LVS_REPORT | LVS_SMALLICON, 0);
	m_cListCtrl.ModifyStyle(0, LVS_REPORT);

	int nColInterval = rect.Width() / 27;
	m_cListCtrl.DeleteColumn(8);
	m_cListCtrl.DeleteColumn(7);
	m_cListCtrl.DeleteColumn(6);
	m_cListCtrl.DeleteColumn(5);
	m_cListCtrl.DeleteColumn(4);
	m_cListCtrl.DeleteColumn(3);
	m_cListCtrl.DeleteColumn(2);
	m_cListCtrl.DeleteColumn(1);
	m_cListCtrl.DeleteColumn(0);

	m_cListCtrl.InsertColumn(0, _T("Set"), LVCFMT_LEFT, nColInterval * 9);
	m_cListCtrl.InsertColumn(1, _T("Size"), LVCFMT_LEFT, 1 * nColInterval);
	m_cListCtrl.InsertColumn(2, _T("Count"), LVCFMT_LEFT, 1 * nColInterval);
	m_cListCtrl.InsertColumn(3, _T("z-score"), LVCFMT_LEFT, 2 * nColInterval);
	m_cListCtrl.InsertColumn(4, _T("Adj. z-score"), LVCFMT_LEFT, 2 * nColInterval);
	m_cListCtrl.InsertColumn(5, _T("Adj. p-value"), LVCFMT_LEFT, 2 * nColInterval);
	m_cListCtrl.InsertColumn(6, _T("q-value"), LVCFMT_LEFT, 2 * nColInterval);
	m_cListCtrl.InsertColumn(7, _T("Network"), LVCFMT_CENTER, 2 * nColInterval);
	m_cListCtrl.InsertColumn(8, _T("List of genes"), LVCFMT_LEFT, rect.Width() - 6 * nColInterval);

	// Insert the default dummy items
	InsertSetItems8();
	m_cListCtrl.SetRedraw(TRUE);
}

void CBigGsaDlg::InitTabCtrl()
{
	m_cTabListMode.InsertItem(0, _T("Set	"));
	m_cTabListMode.EnableWindow(0);
}

void CBigGsaDlg::OnBnClickedRadiohg18() {
	m_fGeneMap.EnableWindow(0);

	m_cCombo_Padding.EnableWindow(1);
	m_cCombo_Padding.SetWindowTextW(_T("±20000"));
	m_rad_hgNumber = 0;
}

void CBigGsaDlg::OnBnClickedRadiohg38()
{
	m_fGeneMap.EnableWindow(0);

	m_cCombo_Padding.EnableWindow(1);
	m_cCombo_Padding.SetWindowTextW(_T("±0"));
	m_rad_hgNumber = 2;
}

void CBigGsaDlg::OnBnClickedRadiohg19()
{
	m_fGeneMap.EnableWindow(0);

	m_cCombo_Padding.EnableWindow(1);
	m_cCombo_Padding.SetWindowTextW(_T("±20000"));
	m_rad_hgNumber = 1;
}

void CBigGsaDlg::OnBnClickedRadiohguser()
{
	m_cCombo_Padding.EnableWindow(0);
	m_fGeneMap.EnableWindow(1);
}

void CBigGsaDlg::OnBnClickedRadioMsigdb()
{
	m_fSetInfo.EnableWindow(0);
	m_cbt_ToSymbol.EnableWindow(0);

	m_cCombo_MSigDB.EnableWindow(1);
	m_cCombo_MSigDB.SetWindowTextW(_T("c2.cp.v5.2.symbols.gmt")); //data\\c5.all.v5.0.symbols.gmt"));
	setFile = _T("c2.cp.v5.2.symbols.gmt");
	//setFile = _T("data\\c5.all.v5.0.symbols.gmt");
}

void CBigGsaDlg::OnBnClickedRadioGo()
{
	m_fSetInfo.EnableWindow(0);
	m_cbt_ToSymbol.EnableWindow(0);
	m_cCombo_MSigDB.EnableWindow(0);

	setFile = _T("data\\Gene Ontology");
}

void CBigGsaDlg::OnBnClickedRadioKegg()
{
	m_fSetInfo.EnableWindow(0);
	m_cbt_ToSymbol.EnableWindow(0);
	m_cCombo_MSigDB.EnableWindow(0);

	setFile = _T("data\\KEGG");
}

void CBigGsaDlg::OnBnClickedRadioOther()
{
	m_cCombo_MSigDB.EnableWindow(0);
	m_fSetInfo.EnableWindow(1);
	m_cbt_ToSymbol.EnableWindow(1);
}

void CBigGsaDlg::OnHdnItemclickListResult(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMHEADER phdr = reinterpret_cast<LPNMHEADER>(pNMHDR);
	// TODO: Add your control notification handler code here
	/*CString str;
	str.Format(_T("%d"),phdr->iItem);
	AfxMessageBox(str);*/
	col_Index = phdr->iItem;

	// Adjusted Z-score
	InsertSetItems8();

	*pResult = 0;
}

void CBigGsaDlg::OnBnClickedRadioSnp()
{
	inputType = SNP_INPUT;
	//m_cCombo_Padding.EnableWindow(1);
	m_fInput.SetWindowTextW(_T("data\\DIAGRAM"));
}

void CBigGsaDlg::OnBnClickedRadioGene()
{
	inputType = GENE_INPUT;
	//m_cCombo_Padding.EnableWindow(0);
	m_fInput.SetWindowTextW(_T("data\\DIAGRAMgene"));
}

void CBigGsaDlg::OnBnClickedbtnperform()
{
	if (test)	test->refreshAll();
	else		test = new gsa();

	////////////////////////////////////////////////////////////////////////////////////////
	m_fInput.GetWindowTextW(inputFile);
	CT2CA pszInputFile(inputFile);
	initParam.inputFile = (string)pszInputFile;

	CString padOption;
	CString dbOption("data\\db");

	if (m_fGeneMap.IsWindowEnabled() != 0)
		m_fGeneMap.GetWindowTextW(mapFile);
	else {
		if (m_rad_hgNumber == 0)		dbOption.Append(_T("18"));
		if (m_rad_hgNumber == 1)		dbOption.Append(_T("19"));
		if (m_rad_hgNumber == 2)		dbOption.Append(_T("38"));

		m_cCombo_Padding.GetWindowTextW(padOption);
		if (padOption.IsEmpty())
			dbOption.Append(_T("_20k"));
		else {
			if (padOption.Compare(CString(_T("±20000"))) == 0)		dbOption.Append(_T("_20k"));
			if (padOption.Compare(CString(_T("±10000"))) == 0)		dbOption.Append(_T("_10k"));
			if (padOption.Compare(CString(_T("0"))) == 0)			dbOption.Append(_T("_0k"));
			if (padOption.Compare(CString(_T("Exon only"))) == 0)	dbOption.Append(_T("_exon"));
		}

		mapFile = dbOption;
		initParam.snpGeneMapFile = (string)(CT2CA)dbOption;
		initParam.hgVersion = m_rad_hgNumber;
	}

	if (m_fSetInfo.IsWindowEnabled() != 0) {
		m_fSetInfo.GetWindowTextW(setFile);
		initParam.setFile = (string)(CT2CA)setFile;
	}

	if (m_cbt_ToSymbol.GetState() == BST_CHECKED) {
		b_Convert2Symbol = 1;
		initParam.convert2Symbol = 1;
	}
	else {
		initParam.convert2Symbol = 0;
		b_Convert2Symbol = 0;
	}

	if (m_cCombo_MSigDB.IsWindowEnabled() != 0)
	{
		CString cOption;
		m_cCombo_MSigDB.GetWindowTextW(cOption);

		if (!cOption.IsEmpty())
		{
			if (cOption.Compare(CString(_T("c1.all.v5.2.symbols.gmt"))) == 0)		setFile = _T("data\\c1.all.v5.2.symbols.gmt");
			if (cOption.Compare(CString(_T("c2.all.v5.2.symbols.gmt"))) == 0)		setFile = _T("data\\c2.all.v5.2.symbols.gmt");
			if (cOption.Compare(CString(_T("c5.all.v5.2.symbols.gmt"))) == 0)		setFile = _T("data\\c5.all.v5.2.symbols.gmt");
			if (cOption.Compare(CString(_T("c2.cp.v5.2.symbols.gmt"))) == 0)		setFile = _T("data\\c2.cp.v5.2.symbols.gmt");
		}

		initParam.setFile = (string)(CT2CA)setFile;
	}

	CString minValue, maxValue;
	m_cEdit_MinGeneSize.GetWindowTextW(minValue);
	if (!minValue.IsEmpty())
	{
		CT2CA pszName(minValue);
		string minSetSize(pszName);
		initParam.minSetSize = atoi(minSetSize.c_str());

	}

	m_cEdit_MaxGeneSize.GetWindowTextW(maxValue);
	if (!maxValue.IsEmpty())
	{
		CT2CA pszName(maxValue);
		string maxSetSize(pszName);
		initParam.maxSetSize = atoi(maxSetSize.c_str());
	}

	//////////////////////////////////////////////////////////////////////////
	CString txt_genecutoff;
	m_Cmb_GeneCutoff.GetWindowTextW(txt_genecutoff);
	if (!txt_genecutoff.IsEmpty())
	{
		CT2CA pszName(txt_genecutoff);
		string genecutoff(pszName);
		this->m_genecutoff = atof(genecutoff.c_str());
		if (this->m_genecutoff <= 0) {
			txt_genecutoff = "0.05";
			this->m_genecutoff = .05;
			initParam.geneCutoff = 0.05;
		}
	}
	else {
		txt_genecutoff = "0.05";
		m_genecutoff = .05;
		initParam.geneCutoff = 0.05;
	}

	CString txt_correlation;
	m_Cmb_AdjacentGeneFile.GetWindowTextW(txt_correlation);
	if (!txt_correlation.IsEmpty())
	{
		if (txt_correlation.Compare(CString(_T("ALL races"))) == 0)			m_correlation = _T("data\\ALL_Adjacent_correlation");
		if (txt_correlation.Compare(CString(_T("African"))) == 0)			m_correlation = _T("data\\AFR_Adjacent_correlation");
		if (txt_correlation.Compare(CString(_T("American"))) == 0)			m_correlation = _T("data\\AMR_Adjacent_correlation");
		if (txt_correlation.Compare(CString(_T("East Asian"))) == 0)		m_correlation = _T("data\\EAS_Adjacent_correlation");
		if (txt_correlation.Compare(CString(_T("European"))) == 0)			m_correlation = _T("data\\EUR_Adjacent_correlation");
		if (txt_correlation.Compare(CString(_T("South Asian"))) == 0)		m_correlation = _T("data\\SAS_Adjacent_correlation");
	}
	else {
		txt_correlation = "European";
		m_correlation = "data\\EUR_Adjacent_correlation";
	}
	initParam.adjacentGeneFile = (string)(CT2CA)m_correlation;

	CString txt_NetworkData;
	m_cbrowser_NetFile.GetWindowTextW(txt_NetworkData);
	if (txt_NetworkData.IsEmpty())
	{
		txt_NetworkData = "data\\STRING_NETWORK.txt";
		m_networkdata = "data\\STRING_NETWORK.txt";		
	}
	else {
		//if (txt_NetworkData.Compare(CString(_T("HIPPIE_NETWORK"))) == 0)			m_networkdata = "data\\HIPPIE_NETWORK.txt";
		//if (txt_NetworkData.Compare(CString(_T("STRING_NETWORK"))) == 0)			m_networkdata = "data\\STRING_NETWORK.txt";
		m_networkdata = txt_NetworkData;
	}
	initParam.networkFile = (string)(CT2CA)m_networkdata;
	initParam.inputType = inputType;

	////////////////////////////////////////////////////////////////////////////////////////
	m_Log.AppendString(_T("Input selected... "));

	if (inputType == SNP_INPUT) {
		m_Log.AppendString(_T("Snp data: ") + inputFile);
		m_Log.AppendString(_T("Gene map data: ") + mapFile);
		m_Log.AppendString(_T("Gene race: ") + txt_correlation);
		m_Log.AppendString(_T("Pathway data: ") + setFile);
		m_Log.AppendString(_T("Pathway size: ") + minValue + _T(" ~ ") + maxValue);
		m_Log.AppendString(_T("Padding: ") + padOption);
		m_Log.AppendString(_T("Additional options..."));
		m_Log.AppendString(_T("Network data: ") + m_networkdata);
		m_Log.AppendString(_T("Network visualization: Gene score cut-off: ") + txt_genecutoff);
		m_Log.AppendString(_T("										"));

		if (inputFile.IsEmpty() ||
			mapFile.IsEmpty() ||
			setFile.IsEmpty())
		{
			m_Log.AppendString(_T("Invalid input file(s)!"));
			m_Log.AppendString(_T("						 "));
			AfxMessageBox(_T("Invalid input file(s)!"));
		}
		else
		{
			m_cListCtrl.EnableWindow(0);
			m_cTabListMode.EnableWindow(0);
			m_fInput.EnableWindow(0);
			m_fGeneMap.EnableWindow(0);
			m_fSetInfo.EnableWindow(0);
			m_cEdit_MaxGeneSize.EnableWindow(0);
			m_cEdit_MinGeneSize.EnableWindow(0);
			m_btnPerform.EnableWindow(0);
			m_EBC_OutputFile.EnableWindow(0);
			m_cCombo_Padding.EnableWindow(0);
			m_cCombo_MSigDB.EnableWindow(0);
			m_cbt_ToSymbol.SetCheck(0);
			m_Cmb_GeneCutoff.EnableWindow(0);
			m_Cmb_AdjacentGeneFile.EnableWindow(0);
			m_CoreNet.EnableWindow(0);
			m_cmb_Net_qvalue_cutoff.EnableWindow(0);
			m_cmb_Net_GeneScoreCutoff.EnableWindow(0);
			m_cbrowser_NetFile.EnableWindow(0);
			m_btnPerform.SetWindowTextW(_T("R U N N I N G"));
			b_Convert2Symbol = 0;

			setWrite = 0;

			/*unique_ptr<mythreadinitparam>	voidParam(make_unique<mythreadinitparam>());
			voidParam->initParam = &initParam;
			voidParam->thisDlg = this;*/
			
			CWinThread *pThread = new CWinThread;			
			pThread = AfxBeginThread(Init8, (PVOID)this, THREAD_PRIORITY_HIGHEST);
			pThread->m_bAutoDelete = FALSE;
			pThread->ResumeThread();
			pThread->Delete();
		}
	}
	else if (inputType == GENE_INPUT) {
		m_Log.AppendString(_T("Gene data: ") + inputFile);
		m_Log.AppendString(_T("Gene map data: ") + mapFile);
		m_Log.AppendString(_T("Gene race: ") + txt_correlation);
		m_Log.AppendString(_T("Pathway data: ") + setFile);
		m_Log.AppendString(_T("Pathway size: ") + minValue + _T(" ~ ") + maxValue);
		m_Log.AppendString(_T("Padding: ") + padOption);
		m_Log.AppendString(_T("Additional options..."));
		m_Log.AppendString(_T("Network data: ") + m_networkdata);
		m_Log.AppendString(_T("Network visualization: Gene score cut-off: ") + txt_genecutoff);
		m_Log.AppendString(_T("										"));

		if (inputFile.IsEmpty() ||
			mapFile.IsEmpty() ||
			setFile.IsEmpty())
		{
			m_Log.AppendString(_T("Invalid input file(s)!"));
			m_Log.AppendString(_T("						 "));
			AfxMessageBox(_T("Invalid input file(s)!"));
		}
		else
		{
			m_cListCtrl.EnableWindow(0);
			m_cTabListMode.EnableWindow(0);
			m_fInput.EnableWindow(0);
			m_fGeneMap.EnableWindow(0);
			m_fSetInfo.EnableWindow(0);
			m_cEdit_MaxGeneSize.EnableWindow(0);
			m_cEdit_MinGeneSize.EnableWindow(0);
			m_btnPerform.EnableWindow(0);
			m_EBC_OutputFile.EnableWindow(0);
			m_cCombo_Padding.EnableWindow(0);
			m_cCombo_MSigDB.EnableWindow(0);
			m_cbt_ToSymbol.SetCheck(0);
			m_Cmb_GeneCutoff.EnableWindow(0);
			m_Cmb_AdjacentGeneFile.EnableWindow(0);
			m_CoreNet.EnableWindow(0);
			m_cmb_Net_qvalue_cutoff.EnableWindow(0);
			m_cmb_Net_GeneScoreCutoff.EnableWindow(0);
			m_cbrowser_NetFile.EnableWindow(0);
			m_btnPerform.SetWindowTextW(_T("R U N N I N G"));
			b_Convert2Symbol = 0;
			setWrite = 0;

			CWinThread *pThread = new CWinThread;
			pThread = AfxBeginThread(Init9, (PVOID)this, THREAD_PRIORITY_HIGHEST);
			pThread->m_bAutoDelete = FALSE;
			pThread->ResumeThread();
			pThread->Delete();
		}
	}
}

void CBigGsaDlg::OnEnUpdateMfcbOutputfile()
{
	// TODO:  If this is a RICHEDIT control, the control will not
	// send this notification unless you override the CDialogEx::OnInitDialog()
	// function to send the EM_SETEVENTMASK message to the control
	// with the ENM_UPDATE flag ORed into the lParam mask.

	// TODO:  Add your control notification handler code here

	UpdateData();
}

void CBigGsaDlg::OnLvnItemchangedListResult(NMHDR *pNMHDR, LRESULT *pResult)
{
	auto_ptr<netparam> param(new netparam);
	param->adjGeneScore = &test->adjGeneScore;
	param->fdrAdjusted = &test->fdrAdjusted;
	param->genenet = &test->genenet;
	param->geneRefinedIndex = &test->geneRefinedIndex;
	param->geneSnpGlobal = &test->geneSnpGlobal;
	param->lgRsInput = &test->lgRsInput;
	param->lgRsInputIndex = &test->lgRsInputIndex;
	param->setGeneRefined = &test->setGeneRefined;
	param->geneInputIndex = &test->geneInputIndex;

	POSITION pos = m_cListCtrl.GetFirstSelectedItemPosition();
	if (pos != NULL) {
		while (pos) {
			int row = m_cListCtrl.GetNextSelectedItem(pos);
			CString setName = m_cListCtrl.GetItemText(row, 0); // Extract set/pathway name
			CT2CA pszName(setName);
			string strSetName(pszName);
			string netfile = test->setSubNetDraw(param.get(), strSetName, m_genecutoff);
			ShellExecute(0, NULL, (LPCTSTR)CString(netfile.c_str()), NULL, NULL, SW_SHOWDEFAULT);
		}
	}
	*pResult = 0;
}

void CBigGsaDlg::OnNMCustomdrawListResult(NMHDR *pNMHDR, LRESULT *pResult)
{
	LPNMLVCUSTOMDRAW lpLVCustomDraw = reinterpret_cast<LPNMLVCUSTOMDRAW>(pNMHDR);

	switch (lpLVCustomDraw->nmcd.dwDrawStage)
	{
	case CDDS_ITEMPREPAINT:
	case CDDS_ITEMPREPAINT | CDDS_SUBITEM:
		//if (0 == ((lpLVCustomDraw->nmcd.dwItemSpec + lpLVCustomDraw->iSubItem) % 2))
		if (lpLVCustomDraw->iSubItem == 7)
		{
			lpLVCustomDraw->clrText = RGB(255, 255, 255); // white text
			lpLVCustomDraw->clrTextBk = RGB(63, 135, 207); // black background
		}
		else
		{
			lpLVCustomDraw->clrText = CLR_DEFAULT;
			lpLVCustomDraw->clrTextBk = CLR_DEFAULT;
		}
		break;

	default: break;
	}

	*pResult = 0;
	*pResult |= CDRF_NOTIFYPOSTPAINT;
	*pResult |= CDRF_NOTIFYITEMDRAW;
	*pResult |= CDRF_NOTIFYSUBITEMDRAW;
}

void CBigGsaDlg::OnBnClickedCorenet()
{
	m_CoreNet.EnableWindow(0);
	m_CoreNet.SetWindowTextW(_T("Drawing"));
	///	
	double netQvalCutoff = 0.25;
	double netGeneCutoff = 0.001;
	//////////////////////////////////////////////////////////////////////////
	CString txt_net_qvalcutoff;
	m_cmb_Net_qvalue_cutoff.GetWindowTextW(txt_net_qvalcutoff);
	if (!txt_net_qvalcutoff.IsEmpty())
	{
		CT2CA pszName(txt_net_qvalcutoff);
		string netqval(pszName);
		double tmp = atof(netqval.c_str());
		if (tmp > 0) {
			netQvalCutoff = tmp;
		}
	}
	else {
		txt_net_qvalcutoff = "0.10";
		netQvalCutoff = 0.10;
	}

	//////////////////////////////////////////////////////////////////////////
	CString txt_net_genecutoff;
	m_cmb_Net_GeneScoreCutoff.GetWindowTextW(txt_net_genecutoff);
	if (!txt_net_qvalcutoff.IsEmpty())
	{
		CT2CA pszName(txt_net_genecutoff);
		string netgene(pszName);
		double tmp = atof(netgene.c_str());
		if (tmp > 0) {
			netGeneCutoff = tmp;
		}
	}
	else {
		txt_net_genecutoff = "0.001";
		netGeneCutoff = 0.001;
	}

	//////////////////////////////////////////////////////////////////////////
	CString outputNetfile;
	m_EBC_OutputFile.GetWindowText(outputNetfile);
	CT2CA pszName(outputNetfile);
	string strNetName(pszName);
	strNetName.insert(0, "core_net_");
	auto_ptr<netparam> param(new netparam);
	param->adjGeneScore = &test->adjGeneScore;
	param->fdrAdjusted = &test->fdrAdjusted;
	param->genenet = &test->genenet;
	param->geneRefinedIndex = &test->geneRefinedIndex;
	param->geneSnpGlobal = &test->geneSnpGlobal;
	param->lgRsInput = &test->lgRsInput;
	param->lgRsInputIndex = &test->lgRsInputIndex;
	param->setGeneRefined = &test->setGeneRefined;
	param->geneInputIndex = &test->geneInputIndex;

	string netfile = test->coreNetDraw(param.get(), strNetName, netQvalCutoff, netGeneCutoff);
	ShellExecute(0, NULL, (LPCTSTR)CString(netfile.c_str()), NULL, NULL, SW_SHOWDEFAULT);

	//////////////////////////////////////////////////////////////////////////
	m_CoreNet.EnableWindow(1);
	m_CoreNet.SetWindowTextW(_T("Global Net"));
}