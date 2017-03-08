
// biggsaDlg.h : header file
//

#pragma once
#include "afxcmn.h"
#include "HistoryEdit.h"
#include "gsa.h"
#include "afxeditbrowsectrl.h"
#include "afxwin.h"


// CBigGsaDlg dialog
class CBigGsaDlg : public CDialogEx
{
// Construction
public:
	CBigGsaDlg(CWnd* pParent = NULL);	// standard constructor
	virtual ~CBigGsaDlg() {
		if (test)	delete test; 
	};
	
// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_BIGGSA_DIALOG };
#endif
		protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support
	
// Implementation
protected:
	HICON m_hIcon;
	
	int col_Index;
	int method_Index;

	// Generated message map functions
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:

	gsa* test;
	CListCtrl	m_cListCtrl;
	CTabCtrl	m_cTabListMode;	
	CHistoryEdit m_Log;
	CMFCEditBrowseCtrl m_fSetInfo;
	CMFCEditBrowseCtrl m_fGeneMap;
	CMFCEditBrowseCtrl m_fInput;
	CButton m_cRadio_hg18;	
	CComboBox m_cCombo_Padding;
	CButton m_cRadio_GO;
	CComboBox m_cCombo_MSigDB;
	CButton m_cRadio_MSigDB;
	CEdit m_cEdit_MinGeneSize;
	CEdit m_cEdit_MaxGeneSize;		
	bool b_Convert2Symbol;

	afx_msg void InitTabCtrl();	
	afx_msg	void InitListCtrlSetView8();
	afx_msg void InsertSetItems8();	
	afx_msg void OnBnClickedRadiohg18();
	afx_msg void OnBnClickedRadiohguser();	
	afx_msg void OnBnClickedRadioGo();
	afx_msg void OnBnClickedRadioOther();	
	afx_msg void OnBnClickedRadioKegg();	
	afx_msg void OnBnClickedRadioMsigdb();
	afx_msg void OnBnClickedRadiohg38();
	afx_msg void OnBnClickedRadiohg19();	
	afx_msg void OnHdnItemclickListResult(NMHDR *pNMHDR, LRESULT *pResult);	
	CButton m_btnPerform;
	afx_msg void OnBnClickedbtnperform();		
	int m_rad_hgNumber;			
	double m_truncation_ratio;	
	double m_genecutoff;	
	CString m_correlation;
	CString m_networkdata;
	int m_NumRS_count;	
	afx_msg void OnEnUpdateMfcbOutputfile();
	CMFCEditBrowseCtrl m_EBC_OutputFile;
	CStatic m_fixed_NumRS;
	afx_msg void OnLvnItemchangedListResult(NMHDR *pNMHDR, LRESULT *pResult);	
	CButton m_cbt_ToSymbol;
	CComboBox m_Cmb_GeneCutoff;
	CComboBox m_Cmb_AdjacentGeneFile;	
	afx_msg void OnNMCustomdrawListResult(NMHDR *pNMHDR, LRESULT *pResult);
	CButton m_CoreNet;
	CComboBox m_cmb_Net_qvalue_cutoff;
	CComboBox m_cmb_Net_GeneScoreCutoff;
	afx_msg void OnBnClickedCorenet();
	CButton m_rad_Input_SNP;
	CButton m_rad_Input_GENE;
	afx_msg void OnBnClickedRadioSnp();
	afx_msg void OnBnClickedRadioGene();
	CMFCEditBrowseCtrl m_cbrowser_NetFile;

};

